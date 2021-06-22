import numpy as np
import scipy.linalg as la

########################################################################################
# Please make use of "./numerical_methods_notes.pdf" to help with understanding the code
########################################################################################

"""
Code Influenced by approaches in:  
1. T. Baumgarte and S. Shapiro, Numerical Relativity. Cambridge: Cambridge University Press, 2010.
2. M. Shibata, Numerical Relativity. Singapore: World Scientific Publishing Co. Pte. Ltd., 2016.
3. T. Baumgarte and S. Shapiro, Numerical Relativity : Starting from Scratch. Cambridge: Cambridge, 2021.
"""


class EllipticEquationSolver:
    """This class solves Poisson-type elliptic equations of the form:
          D^2 solution + fct solution = rhs
    where
        - D^2 is the flat Laplace operator
        - fct and rhs are user-supplied functions of the coordinates x, y, z,
        - and solution is the solution.
    """

    def __init__(self, x: float, y: float, z: float) -> None:
        """Constructor - provide Cartesian coordinates, all of length n_grid,
        as arguments.
        """

        print("Constructing Poisson solver...")
        self.n_grid = np.size(x)
        self.delta = x[1] - x[0]

        # set up storage for matrix, solution, r.h.s.
        # Note: "solution" and "RHS" will store functions in 3d format, while
        # "solution_1D" and "RHS_1D" will store functions in 1d format
        nnn = self.n_grid ** 3
        self.RHS_1D = np.zeros(nnn)
        self.A = np.zeros((nnn, nnn))
        self.solution = np.zeros((self.n_grid, self.n_grid, self.n_grid))
        self.radius = np.zeros((self.n_grid, self.n_grid, self.n_grid))

        # compute radius
        for i in range(0, self.n_grid):
            for j in range(0, self.n_grid):
                for k in range(0, self.n_grid):
                    r2 = x[i] ** 2 + y[j] ** 2 + z[k] ** 2
                    self.radius[i, j, k] = np.sqrt(r2)

    def setup_matrix(self, fct: np.ndarray) -> None:
        """Set up matrix A."""

        n_grid = self.n_grid

        # Use Robin boundary conditions to set up boundaries
        i = 0  # lower x-boundary
        for j in range(0, n_grid):
            for k in range(0, n_grid):
                super_index = self.__super_index(i, j, k)
                self.A[super_index, super_index] = self.radius[i, j, k]
                self.A[super_index, super_index + 1] = - \
                    self.radius[i + 1, j, k]

        i = n_grid - 1  # upper x-boundary
        for j in range(0, n_grid):
            for k in range(0, n_grid):
                super_index = self.__super_index(i, j, k)
                self.A[super_index, super_index] = self.radius[i, j, k]
                self.A[super_index, super_index - 1] = - \
                    self.radius[i - 1, j, k]

        j = 0  # lower y-boundary
        for i in range(1, n_grid - 1):
            for k in range(0, n_grid):
                super_index = self.__super_index(i, j, k)
                self.A[super_index, super_index] = self.radius[i, j, k]
                self.A[super_index, super_index + n_grid] = - \
                    self.radius[i, j + 1, k]

        j = n_grid - 1  # upper y-boundary
        for i in range(1, n_grid - 1):
            for k in range(0, n_grid):
                super_index = self.__super_index(i, j, k)
                self.A[super_index, super_index] = self.radius[i, j, k]
                self.A[super_index, super_index - n_grid] = - \
                    self.radius[i, j - 1, k]

        k = 0  # lower z-boundary
        for i in range(1, n_grid - 1):
            for j in range(1, n_grid - 1):
                super_index = self.__super_index(i, j, k)
                self.A[super_index, super_index] = self.radius[i, j, k]
                self.A[super_index, super_index + n_grid * n_grid] = - \
                    self.radius[i, j, k + 1]

        k = n_grid - 1  # upper z-boundary
        for i in range(1, n_grid - 1):
            for j in range(1, n_grid - 1):
                super_index = self.__super_index(i, j, k)
                self.A[super_index, super_index] = self.radius[i, j, k]
                self.A[super_index, super_index - n_grid * n_grid] = - \
                    self.radius[i, j, k - 1]

        # fill matrix in interior
        for i in range(1, n_grid - 1):
            for j in range(1, n_grid - 1):
                for k in range(1, n_grid - 1):
                    super_index = self.__super_index(i, j, k)

                    # diagonal element
                    self.A[super_index, super_index] = - \
                        6. + self.delta ** 2 * fct[i, j, k]

                    # off-diagonal elements
                    self.A[super_index, super_index - 1] = 1.0
                    self.A[super_index, super_index + 1] = 1.0
                    self.A[super_index, super_index - n_grid] = 1.0
                    self.A[super_index, super_index + n_grid] = 1.0
                    self.A[super_index, super_index - n_grid * n_grid] = 1.0
                    self.A[super_index, super_index + n_grid * n_grid] = 1.0

    def setup_rhs(self, RHS: np.ndarray) -> None:
        """Setup right-hand side of matrix equation"""

        n_grid = self.n_grid
        for i in range(1, n_grid - 1):
            for j in range(1, n_grid - 1):
                for k in range(1, n_grid - 1):
                    super_index = self.__super_index(i, j, k)
                    self.RHS_1D[super_index] = self.delta ** 2 * RHS[i, j, k]

    def solve(self) -> np.ndarray:
        """Interface to scipy.linalg matrix solver,
        returns solution (in 3d format)."""

        # solve matrix using scipy.linalg interface...
        solution_1D = la.solve(self.A, self.RHS_1D)

        # ... then translate from superindex to 3d
        for i in range(0, self.n_grid):
            for j in range(0, self.n_grid):
                for k in range(0, self.n_grid):
                    super_index = self.__super_index(i, j, k)
                    self.solution[i, j, k] = solution_1D[super_index]

        return self.solution

    def __super_index(self, i: int, j: int, k: int) -> int:
        """Compute super index"""
        return i + self.n_grid * (j + self.n_grid * k)
