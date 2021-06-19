import numpy as np
import scipy.linalg as la

########################################################################################
# Please make use of "./numerical_methods_notes.pdf" to help with understanding the code

# Code reference : T. Baumgarte and S. Shapiro, Numerical relativity. Cambridge: Cambridge University Press, 2010.
########################################################################################


class EllipticSolver:
    """Class Elliptic solves Poisson-type elliptic equations of the form:
          D^2 sol + fct sol = rhs
    where
        - D^2 is the flat Laplace operator
        - fct and rhs are user-supplied functions of the coordinates x, y, z,
        - and sol is the solution.
    """

    def __init__(self, x: float, y: float, z: float) -> None:
        """Constructor - provide Cartesian coordinates, all of length n_grid,
        as arguments.
        """

        print("Setting up Poisson solver...")
        self.n_grid = np.size(x)
        self.delta = x[1] - x[0]

        # set up storage for matrix, solution, r.h.s.
        # Note: "sol" and "rhs" will store functions in 3d format, while
        # "sol_1d" and "rhs_1d" will store functions in 1d format
        nnn = self.n_grid ** 3
        self.rhs_1d = np.zeros(nnn)
        self.A = np.zeros((nnn, nnn))
        self.sol = np.zeros((self.n_grid, self.n_grid, self.n_grid))
        self.rad = np.zeros((self.n_grid, self.n_grid, self.n_grid))

        # compute radius
        for i in range(0, self.n_grid):
            for j in range(0, self.n_grid):
                for k in range(0, self.n_grid):
                    rad2 = x[i] ** 2 + y[j] ** 2 + z[k] ** 2
                    self.rad[i, j, k] = np.sqrt(rad2)

    def setup_matrix(self, fct: np.ndarray) -> None:
        """Set up matrix A."""

        n_grid = self.n_grid

        # Use Robin boundary conditions to set up boundaries
        i = 0  # lower x-boundary
        for j in range(0, n_grid):
            for k in range(0, n_grid):
                index = self.__super_index(i, j, k)
                self.A[index, index] = self.rad[i, j, k]
                self.A[index, index + 1] = -self.rad[i + 1, j, k]

        i = n_grid - 1  # upper x-boundary
        for j in range(0, n_grid):
            for k in range(0, n_grid):
                index = self.__super_index(i, j, k)
                self.A[index, index] = self.rad[i, j, k]
                self.A[index, index - 1] = -self.rad[i - 1, j, k]

        j = 0  # lower y-boundary
        for i in range(1, n_grid - 1):
            for k in range(0, n_grid):
                index = self.__super_index(i, j, k)
                self.A[index, index] = self.rad[i, j, k]
                self.A[index, index + n_grid] = -self.rad[i, j + 1, k]

        j = n_grid - 1  # upper y-boundary
        for i in range(1, n_grid - 1):
            for k in range(0, n_grid):
                index = self.__super_index(i, j, k)
                self.A[index, index] = self.rad[i, j, k]
                self.A[index, index - n_grid] = -self.rad[i, j - 1, k]

        k = 0  # lower z-boundary
        for i in range(1, n_grid - 1):
            for j in range(1, n_grid - 1):
                index = self.__super_index(i, j, k)
                self.A[index, index] = self.rad[i, j, k]
                self.A[index, index + n_grid * n_grid] = -self.rad[i, j, k + 1]

        k = n_grid - 1  # upper z-boundary
        for i in range(1, n_grid - 1):
            for j in range(1, n_grid - 1):
                index = self.__super_index(i, j, k)
                self.A[index, index] = self.rad[i, j, k]
                self.A[index, index - n_grid * n_grid] = -self.rad[i, j, k - 1]

        # fill matrix in interior
        for i in range(1, n_grid - 1):
            for j in range(1, n_grid - 1):
                for k in range(1, n_grid - 1):
                    index = self.__super_index(i, j, k)

                    # diagonal element
                    self.A[index, index] = -6. + self.delta ** 2 * fct[i, j, k]

                    # off-diagonal elements
                    self.A[index, index - 1] = 1.0
                    self.A[index, index + 1] = 1.0
                    self.A[index, index - n_grid] = 1.0
                    self.A[index, index + n_grid] = 1.0
                    self.A[index, index - n_grid * n_grid] = 1.0
                    self.A[index, index + n_grid * n_grid] = 1.0

    def setup_rhs(self, rhs: np.ndarray) -> None:
        """Setup right-hand side of matrix equation"""

        n_grid = self.n_grid
        for i in range(1, n_grid - 1):
            for j in range(1, n_grid - 1):
                for k in range(1, n_grid - 1):
                    index = self.__super_index(i, j, k)
                    self.rhs_1d[index] = self.delta ** 2 * rhs[i, j, k]

    def solve(self) -> np.ndarray:
        """Interface to scipy.linalg matrix solver,
        returns sol (in 3d format)."""

        # solve matrix using scipy.linalg interface...
        sol_1d = la.solve(self.A, self.rhs_1d)

        # ... then translate from superindex to 3d
        for i in range(0, self.n_grid):
            for j in range(0, self.n_grid):
                for k in range(0, self.n_grid):
                    index = self.__super_index(i, j, k)
                    self.sol[i, j, k] = sol_1d[index]

        return self.sol

    def __super_index(self, i: int, j: int, k: int) -> int:
        """Compute super index"""
        return i + self.n_grid * (j + self.n_grid * k)
