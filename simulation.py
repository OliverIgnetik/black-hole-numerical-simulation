"""Code to construct puncture initial data for single black hole."""
import sys
from numpy import zeros, sqrt, linspace

from elliptic_solver import EllipticSolver

class Puncture:
    """Class that handles construction of puncture data.

    To use this class,
       - initialize class with physical parameters as arguments
       - then call construct_solution.
    """

    def __init__(self, bh_loc, lin_mom, n_grid, x_out):
        """Arguments to constructor specify physical parameters:
        - location of puncture (bh_loc)
        - linear momentum (lin_mom)
        - size of grid (n_grid)
        - outer boundary (x_out).
        """
        self.bh_loc = bh_loc
        self.lin_mom = lin_mom
        # echo out parameters
        print(" Constructing class Puncture for single black hole")
        print("    at bh_loc = (", bh_loc[0], ",", bh_loc[1], ",",
              bh_loc[2], ")")
        print("    with momentum p = (", lin_mom[0], ",",
              lin_mom[1], ",", lin_mom[2], ")")
        print(" Using", n_grid, "\b^3 gridpoints with outer boundary at", x_out)
        # set up grid
        self.n_grid = n_grid
        self.x_out = x_out
        self.delta = 2.0 * x_out / n_grid

        # set up coordinates: use cell-centered grid covering (-x_out, x_out)
        # in each dimension; see (B.14)
        half_delta = self.delta / 2.0
        self.x = linspace(half_delta - x_out, x_out -
                          half_delta, n_grid)
        self.y = linspace(half_delta - x_out, x_out -
                          half_delta, n_grid)
        self.z = linspace(half_delta - x_out, x_out -
                          half_delta, n_grid)

        # allocate elliptic solver
        self.solver = EllipticSolver(self.x, self.y, self.z)

        # allocate memory for functions u, alpha, beta, and residual
        self.alpha = zeros((n_grid, n_grid, n_grid))
        self.beta = zeros((n_grid, n_grid, n_grid))
        self.u = zeros((n_grid, n_grid, n_grid))
        self.res = zeros((n_grid, n_grid, n_grid))

    def construct_solution(self, tol, it_max):
        """Construct solution iteratively, provide tolerance and maximum
        number of iterations as arguments."""

        self.setup_alpha_beta()
        residual_norm = self.residual()
        print(" Initial Residual = ", residual_norm)
        print(" Using up to", it_max, "iteration steps to reach tolerance of",
              tol)

        # now iterate...
        it_step = 0
        while residual_norm > tol and it_step < it_max:
            it_step += 1
            self.update_u()
            residual_norm = self.residual()
            print(" Residual after", it_step, "iterations :", residual_norm)
        if (residual_norm < tol):
            print(" Done!")
        else:
            print(" Giving up...")

    def update_u(self):
        """Function that updates u using Poisson solver;
        takes one iteration step.
        """

        # set up linear term and right-hand side for SolvePoisson...
        n_grid = self.n_grid
        fct = zeros((n_grid, n_grid, n_grid))
        rhs = zeros((n_grid, n_grid, n_grid))

        for i in range(1, n_grid - 1):
            for j in range(1, n_grid - 1):
                for k in range(1, n_grid - 1):
                    # compute h' from (B.39)
                    temp = self.alpha[i, j, k] * (1.0 + self.u[i, j, k]) + 1.0
                    fct[i, j, k] = (-7.0 * self.beta[i, j, k] *
                                    self.alpha[i, j, k] / temp ** 8)
                    rhs[i, j, k] = -self.res[i, j, k]

        # now update Poisson solver
        self.solver.setup_matrix(fct)

        # set up right-hand side
        self.solver.setup_rhs(rhs)

        # solve to find delta_u, see (B.36)
        delta_u = self.solver.solve()

        # update u
        self.u += delta_u

    def residual(self):
        """Evaluate residual, see (B.35)."""

        residual_norm = 0.0
        for i in range(1, self.n_grid - 1):
            for j in range(1, self.n_grid - 1):
                for k in range(1, self.n_grid - 1):

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

                    # then compute difference to get residual, see (B.35)
                    self.res[i, j, k] = lhs - rhs
                    residual_norm += self.res[i, j, k] ** 2

        residual_norm = sqrt(residual_norm) * self.delta ** 3
        return residual_norm

    def setup_alpha_beta(self):
        """Set up functions alpha and beta."""

        n_grid = self.n_grid
        p_x = self.lin_mom[0]
        p_y = self.lin_mom[1]
        p_z = self.lin_mom[2]

        for i in range(0, n_grid):
            for j in range(0, n_grid):
                for k in range(0, n_grid):
                    s_x = self.x[i] - self.bh_loc[0]
                    s_y = self.y[j] - self.bh_loc[1]
                    s_z = self.z[k] - self.bh_loc[2]
                    s2 = s_x ** 2 + s_y ** 2 + s_z ** 2
                    s_bh = sqrt(s2)
                    l_x = s_x / s_bh
                    l_y = s_y / s_bh
                    l_z = s_z / s_bh
                    lP = l_x * p_x + l_y * p_y + l_z * p_z

                    # construct extrinsic curvature, see Eq. (3.43)
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

                    # now compute alpha and beta from (3.47) and (3.49)
                    self.alpha[i, j, k] = 2.0 * s_bh
                    self.beta[i, j, k] = self.alpha[i, j, k] ** 7 * A2 / 8.0

    def write_to_file(self):
        """Function that writes solution to file."""

        filename = "simulation_data_" + str(self.n_grid) + "_" + str(self.x_out)
        filename = filename + ".data"
        out = open(filename, "w")
        if out:
            k = self.n_grid // 2
            out.write(
                "# Data for black hole at x = (%f,%f,%f)\n"
                % (self.bh_loc[0], self.bh_loc[1], self.bh_loc[2])
            )
            out.write("# with linear momentum P = (%f, %f, %f)\n" %
                      (self.lin_mom))
            out.write("# in plane for z = %e \n" % (self.z[k]))
            out.write("# x            y              u              \n")
            out.write("#============================================\n")
            for i in range(0, self.n_grid):
                for j in range(0, self.n_grid):
                    out.write("%e  %e  %e\n" % (self.x[i], self.y[j],
                                                self.u[i, j, k]))
                out.write("\n")
            out.close()
        else:
            print(" Could not open file", filename, "in write_to_file()")
            print(" Check permissions?")
#
# =====================================================================
# Main routine: defines parameters, sets up puncture solver, and
# then finds solution
# =====================================================================
#


def main():
    """Main routine..."""
    print(" -------------------------------------------------------")
    print(" --- puncture.py --- use flag -h for list of options ---")
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
    n_grid = 16
    # location of outer boundary
    x_out = 4.0
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
        if sys.argv[i] == "-n_grid":
            n_grid = int(sys.argv[i+1])
        if sys.argv[i] == "-x_out":
            x_out = float(sys.argv[i+1])
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
    bh_loc = (loc_x, loc_y, loc_z)
    # linear momentum
    lin_mom = (p_x, p_y, p_z)
    #
    # set up Puncture solver
    black_hole = Puncture(bh_loc, lin_mom, n_grid, x_out)
    #
    # and construct solution
    black_hole.construct_solution(tol, it_max)
    #
    # and write results to file
    black_hole.write_to_file()


def usage():
    print("Constructs puncture initial data for single black hole.")
    print("")
    print("The following options can be used to over-write default parameters")
    print("\t-n_grid: number of grid points [default: 16]")
    print("\t-x_out: location of outer boundary [4.0]")
    print("\t-loc_x, -loc_y, -loc_z: location of black hole [(0.0, 0.0, 0.0)]")
    print("\t-p_x, -p_y, -p_z: lin. momentum of black hole [(1.0, 0.0, 0.0)]")
    print("\t-tol: tolerance for elliptic solver [1.e-12]")
    print("\t-it_max: maximum number of iterations [50]")
    print("For example, to construct data with x_out = 6.0, call")
    print("\tpython3 puncture.py -x_out 6.0")


if __name__ == '__main__':
    main()
