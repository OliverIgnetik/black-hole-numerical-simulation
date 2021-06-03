import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
import numpy as np
import matplotlib
from distutils.spawn import find_executable


def puncture_plot(data_file, plot_file=0):
    if find_executable('latex'):
        rc('text', usetex=True)
    #
    f = open(data_file, 'r')
    if f:
        print("Reading data from file", data_file)
    else:
        print("Cannot open data file", data_file)
        return
    x, y, fct = np.loadtxt(data_file, unpack=True)
    f.close()
    #
    n_grid = int(np.sqrt(x.size))
    X = np.reshape(x, (n_grid, n_grid))
    Y = np.reshape(y, (n_grid, n_grid))
    FCT = np.reshape(fct, (n_grid, n_grid))
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, FCT, cmap=cm.autumn, linewidth=1,
                    antialiased=False, alpha=0.2)
    ax.contour(X, Y, FCT, zdir='z', offset=0.0, cmap=cm.autumn)
    #
    if find_executable('latex'):
        ax.set_xlabel(r"$x / {\mathcal M}$", size=14)
        ax.set_ylabel(r"$y / {\mathcal M}$", size=14)
        ax.set_zlabel(r"$u$", size=14)
    else:
        ax.set_xlabel("x / M", size=14)
        ax.set_ylabel("y / M", size=14)
        ax.set_zlabel("u", size=14)
    ax.set_zlim(0.0, max(fct))
    ax.tick_params(axis='x', which='major', pad=-2)
    ax.tick_params(axis='y', which='major', pad=-2)
    ax.xaxis.set_rotate_label(False)
    ax.yaxis.set_rotate_label(False)
    ax.zaxis.set_rotate_label(False)
    ax.view_init(elev=20., azim=-110)
    if plot_file == 0:
        plt.show()
    else:
        plt.savefig(plot_file, psi=300)


def main():
    """Main routine..."""
    print(" ------------------------------------------------------------")
    print(" --- puncture_plot.py --- use flag -h for list of options ---")
    print(" ------------------------------------------------------------")
    # print("Using matplotlib version", matplotlib.__version__)
    #
    # set default values for variables
    #
    # data filename
    data_file = "Puncture_16_4.0.data"
    # filename to be plotted to (set to zero for display on screen)
    plot_file = 0
    #
    # now look for flags to overwrite default values
    #
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-h":
            usage()
            return
        if sys.argv[i] == "-data":
            data_file = sys.argv[i+1]
        if sys.argv[i] == "-plot":
            plot_file = sys.argv[i+1]
    #
    # plot puncture data...
    #
    puncture_plot(data_file, plot_file)


def usage():
    print("Creates plots from data files produced with puncture.py.")
    print("")
    print("The following options can be used to over-write default parameters")
    print(
        "\t-data: provide name of data_file [Default: Puncture_16_4.0.data] ")
    print(
        "\t-plot: if provided, plot will be saved to this file [Default: None]")
    print("")
    print("For example, to make a plot of the file Puncture_24_6.0.data")
    print("and save to file plot.pdf, call")
    print("\tpython3 puncture_plot.py -data Puncture_24_6.0.data -plot plot.pdf")


if __name__ == '__main__':
    main()


# puncture()
# compare_formulations()
# convergence()
# outbound_convergence()
