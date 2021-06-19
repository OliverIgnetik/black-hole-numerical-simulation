"""Code to construct puncture initial data for single black hole."""
import sys
from numpy import zeros, sqrt, linspace
from elliptic_solver import EllipticEquationSolver
from typing import Tuple
"""
Code References
1. T. Baumgarte and S. Shapiro, Numerical Relativity. Cambridge: Cambridge University Press, 2010.
2. T. Baumgarte and S. Shapiro, Numerical Relativity : Starting from Scratch. Cambridge: Cambridge, 2021.
"""

########################################################################################
# Please make use of "./numerical_methods_notes.pdf" to help with understanding the code
########################################################################################


class Puncture:
    """Class that handles construction of puncture data."""

    def __init__(self, bh_location: Tuple[float, float, float], linear_momentum: Tuple[float, float, float], grid_dim: int, boundary: float) -> None:
        """Arguments to constructor specify physical parameters:
        - location of puncture (bh_location)
        - linear momentum (linear_momentum)
        - size of grid (grid_dim)
        - outer boundary (boundary).
        """
        self.bh_location = bh_location
        self.linear_momentum = linear_momentum
        # echo out parameters
        print("Constructing puncture data for single black hole")
        print(
            f"at bh_location = ({bh_location[0]}, {bh_location[1]}, {bh_location[2]} )")
        print(
            f"with momentum p = ({linear_momentum[0]},{linear_momentum[1]}, {linear_momentum[2]}")
        print(
            f"Using {grid_dim}\b^3 gridpoints with outer boundary at {boundary}")
        # set up grid
        self.grid_dim = grid_dim
        self.boundary = boundary
        self.delta = 2.0 * boundary / grid_dim

        # set up coordinates: use cell-centered grid covering (-boundary, boundary) in each dimension
        half_delta = self.delta / 2.0
        self.x = linspace(half_delta - boundary, boundary -
                          half_delta, grid_dim)
        self.y = linspace(half_delta - boundary, boundary -
                          half_delta, grid_dim)
        self.z = linspace(half_delta - boundary, boundary -
                          half_delta, grid_dim)

        # setup elliptic solver
        self.solver = EllipticEquationSolver(self.x, self.y, self.z)

        # allocate memory for functions u, alpha, beta, and residual
        self.alpha = zeros((grid_dim, grid_dim, grid_dim))
        self.beta = zeros((grid_dim, grid_dim, grid_dim))
        self.u = zeros((grid_dim, grid_dim, grid_dim))
        self.res = zeros((grid_dim, grid_dim, grid_dim))

    def construct_solution(self, tol: float, it_max: float) -> None:
        """Construct solution iteratively, provide tolerance and maximum
        number of iterations as arguments."""

        self.__setup_alpha_beta()
        residual_norm = self.__residual()
        print(f"Initial Residual = { residual_norm}")
        print(
            f" Maximum iterations : {it_max} steps to reach tolerance : {tol}")

        # now iterate...
        it_step = 0
        while residual_norm > tol and it_step < it_max:
            it_step += 1
            self.__update_u()
            residual_norm = self.__residual()
            print(f"Residual : {it_step} iterations : {residual_norm}")

        if (residual_norm < tol):
            print(" Tolerance level for residual achieved.")
        else:
            raise RuntimeError(
                'It is not possible to reduce the residual norm below {tol}\n Terminating simulation...')

    def __update_u(self) -> None:
        """Function that updates u using Poisson solver;
        takes one iteration step.
        """

        # set up linear term and right-hand side for SolvePoisson...
        grid_dim = self.grid_dim
        fct = zeros((grid_dim, grid_dim, grid_dim))
        rhs = zeros((grid_dim, grid_dim, grid_dim))

        for i in range(1, grid_dim - 1):
            for j in range(1, grid_dim - 1):
                for k in range(1, grid_dim - 1):
                    # compute h'
                    temp = self.alpha[i, j, k] * (1.0 + self.u[i, j, k]) + 1.0
                    fct[i, j, k] = (-7.0 * self.beta[i, j, k] *
                                    self.alpha[i, j, k] / temp ** 8)
                    rhs[i, j, k] = -self.res[i, j, k]

        # now update Poisson solver
        self.solver.setup_matrix(fct)

        # set up right-hand side
        self.solver.setup_rhs(rhs)

        # solve to find delta_u
        delta_u = self.solver.solve()

        # update u
        self.u += delta_u

    def __residual(self) -> float:
        """Evaluate residual"""

        residual_norm = 0.0
        for i in range(1, self.grid_dim - 1):
            for j in range(1, self.grid_dim - 1):
                for k in range(1, self.grid_dim - 1):

                    # compute left-hand side: Laplace operator
                    ddx = (self.u[i + 1, j, k] - 2.0 * self.u[i, j, k] +
                           self.u[i - 1, j, k])
                    ddy = (self.u[i, j + 1, k] - 2.0 * self.u[i, j, k] +
                           self.u[i, j - 1, k])
                    ddz = (self.u[i, j, k + 1] - 2.0 * self.u[i, j, k] +
                           self.u[i, j, k - 1])
                    lhs = (ddx + ddy + ddz) / self.delta ** 2

                    # compute right-hand side,
                    # recall h = - beta/(alpha + alpha u + 1)^7
                    temp = self.alpha[i, j, k] * (1.0 + self.u[i, j, k]) + 1.0
                    rhs = -self.beta[i, j, k] / temp ** 7

                    # then compute difference to get residual
                    self.res[i, j, k] = lhs - rhs
                    residual_norm += self.res[i, j, k] ** 2

        residual_norm = sqrt(residual_norm) * self.delta ** 3
        return residual_norm

    def __setup_alpha_beta(self) -> None:
        """Set up functions alpha and beta."""

        grid_dim = self.grid_dim
        p_x = self.linear_momentum[0]
        p_y = self.linear_momentum[1]
        p_z = self.linear_momentum[2]

        for i in range(0, grid_dim):
            for j in range(0, grid_dim):
                for k in range(0, grid_dim):
                    s_x = self.x[i] - self.bh_location[0]
                    s_y = self.y[j] - self.bh_location[1]
                    s_z = self.z[k] - self.bh_location[2]
                    s2 = s_x ** 2 + s_y ** 2 + s_z ** 2
                    s_bh = sqrt(s2)
                    l_x = s_x / s_bh
                    l_y = s_y / s_bh
                    l_z = s_z / s_bh
                    lP = l_x * p_x + l_y * p_y + l_z * p_z

                    # construct extrinsic curvature
                    fac = 3.0 / (2.0 * s2)
                    A_xx = fac * (2.0 * p_x * l_x - (1.0 - l_x * l_x) * lP)
                    A_yy = fac * (2.0 * p_y * l_y - (1.0 - l_y * l_y) * lP)
                    A_zz = fac * (2.0 * p_z * l_z - (1.0 - l_z * l_z) * lP)
                    A_xy = fac * (p_x * l_y + p_y * l_x + l_x * l_y * lP)
                    A_xz = fac * (p_x * l_z + p_z * l_x + l_x * l_z * lP)
                    A_yz = fac * (p_y * l_z + p_z * l_y + l_y * l_z * lP)

                    # compute A_{ij} A^{ij}
                    A2 = (
                        A_xx ** 2 + A_yy ** 2 + A_zz ** 2 +
                        2.0*(A_xy ** 2 + A_xz ** 2 + A_yz ** 2)
                    )

                    # now compute alpha and beta
                    self.alpha[i, j, k] = 2.0 * s_bh
                    self.beta[i, j, k] = self.alpha[i, j, k] ** 7 * A2 / 8.0

    def write_to_file(self) -> None:
        """Function that writes solution to file."""

        filename = "simulation_data_" + \
            str(self.grid_dim) + "_" + str(self.boundary)
        filename = filename + ".data"
        out = open(filename, "w")
        if out:
            k = self.grid_dim // 2
            out.write(
                "# Data for black hole at x = (%f,%f,%f)\n"
                % (self.bh_location[0], self.bh_location[1], self.bh_location[2])
            )
            out.write("# with linear momentum P = (%f, %f, %f)\n" %
                      (self.linear_momentum))
            out.write("# in plane for z = %e \n" % (self.z[k]))
            out.write("# x            y              u              \n")
            out.write("#============================================\n")
            for i in range(0, self.grid_dim):
                for j in range(0, self.grid_dim):
                    out.write("%e  %e  %e\n" % (self.x[i], self.y[j],
                                                self.u[i, j, k]))
                out.write("\n")
            out.close()
        else:
            raise PermissionError(
                f"Could not open file filename, in write_to_file() \nCheck permissions?")
#
# =====================================================================
# Main routine: defines parameters, sets up puncture solver, and
# then finds solution
# =====================================================================
#


def main():
    """Main routine..."""
    print(" -------------------------------------------------------")
    print(" --- simulation.py --- use flag -h for list of options ---")
    print(" -------------------------------------------------------")
    #
    # set default values for variables
    #
    # location of black hole:
    loc_x = 0.0
    loc_y = 0.0
    loc_z = 0.0
    # momentum of black hole:
    p_x = 1.0
    p_y = 0.0
    p_z = 0.0
    # number of grid points
    grid_dim = 16
    # location of outer boundary
    boundary = 4.0
    # tolerance and maximum number of iterations
    tol = 1.0e-12
    it_max = 50
    #
    # now look for flags to overwrite default values
    #
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-h":
            usage()
            return
        if sys.argv[i] == "-grid_dim":
            grid_dim = int(sys.argv[i+1])
        if sys.argv[i] == "-boundary":
            boundary = float(sys.argv[i+1])
        if sys.argv[i] == "-loc_x":
            loc_x = float(sys.argv[i+1])
        if sys.argv[i] == "-loc_y":
            loc_y = float(sys.argv[i+1])
        if sys.argv[i] == "-loc_z":
            loc_z = float(sys.argv[i+1])
        if sys.argv[i] == "-p_x":
            p_x = float(sys.argv[i+1])
        if sys.argv[i] == "-p_y":
            p_y = float(sys.argv[i+1])
        if sys.argv[i] == "-p_z":
            p_z = float(sys.argv[i+1])
        if sys.argv[i] == "-tol":
            tol = float(sys.argv[i+1])
        if sys.argv[i] == "-it_max":
            it_max = int(sys.argv[i+1])

    # location of puncture
    bh_location = (loc_x, loc_y, loc_z)
    # linear momentum
    linear_momentum = (p_x, p_y, p_z)
    #
    # set up Puncture solver
    black_hole = Puncture(bh_location, linear_momentum, grid_dim, boundary)
    #
    # and construct solution
    black_hole.construct_solution(tol, it_max)
    #
    # and write results to file
    black_hole.write_to_file()


def usage():
    print("Constructs puncture initial data for single black hole.")
    print("-"*30)
    print("The following options can be used to over-write default parameters")
    print("\t-grid_dim: number of grid points [default: 16]")
    print("\t-boundary: location of outer boundary [4.0]")
    print("\t-loc_x, -loc_y, -loc_z: location of black hole [(0.0, 0.0, 0.0)]")
    print("\t-p_x, -p_y, -p_z: lin. momentum of black hole [(1.0, 0.0, 0.0)]")
    print("\t-tol: tolerance for elliptic solver [1.e-12]")
    print("\t-it_max: maximum number of iterations [50]")
    print("For example, to construct data with boundary = 6.0, call")
    print("\tpython simulation.py -boundary 6.0")


if __name__ == '__main__':
    main()
