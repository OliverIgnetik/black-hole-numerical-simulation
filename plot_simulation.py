import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
import numpy as np
import matplotlib
from distutils.spawn import find_executable

"""
Code Influenced by approaches in:  
1. T. Baumgarte and S. Shapiro, Numerical Relativity. Cambridge: Cambridge University Press, 2010.
2. M. Shibata, Numerical Relativity. Singapore: World Scientific Publishing Co. Pte. Ltd., 2016.
3. T. Baumgarte and S. Shapiro, Numerical Relativity : Starting from Scratch. Cambridge: Cambridge, 2021.
"""


def puncture_plot(data_file: str, plot_file=0) -> None:
    if find_executable('latex'):
        rc('text', usetex=True)
    #
    f = open(data_file, 'r')
    if f:
        print("Reading data from file", data_file)
    else:
        raise FileNotFoundError(f"Cannot open data file {data_file}")
    x, y, fct = np.loadtxt(data_file, unpack=True)
    f.close()
    #
    n_grid = int(np.sqrt(x.size))
    X = np.reshape(x, (n_grid, n_grid))
    Y = np.reshape(y, (n_grid, n_grid))
    FCT = np.reshape(fct, (n_grid, n_grid))
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, FCT, cmap=cm.plasma, linewidth=1,
                    antialiased=False, alpha=0.2)
    ax.contour(X, Y, FCT, zdir='z', offset=0.0, cmap=cm.plasma)
    #
    if find_executable('latex'):
        ax.set_xlabel(r"$x / {\mathcal M}$", size=14)
        ax.set_ylabel(r"$y / {\mathcal M}$", size=14)
        ax.set_zlabel(r"$U$", size=14)
    else:
        ax.set_xlabel("x / M", size=14)
        ax.set_ylabel("y / M", size=14)
        ax.set_zlabel("U ", size=14)
    ax.set_zlim(0.0, max(fct))
    ax.tick_params(axis='x', which='major', pad=-2)
    ax.tick_params(axis='y', which='major', pad=-2)
    ax.xaxis.set_rotate_label(False)
    ax.yaxis.set_rotate_label(False)
    ax.zaxis.set_rotate_label(False)
    ax.view_init(elev=23., azim=-100)
    if plot_file == 0:
        plt.show()
    else:
        plt.savefig(plot_file)


def main():
    """Main routine..."""
    print(" ------------------------------------------------------------")
    print(" --- simulation_plot.py --- use flag -h for list of options ---")
    print(" ------------------------------------------------------------")
    # print("Using matplotlib version", matplotlib.__version__)
    #
    # set default values for variables
    #
    # data filename
    data_file = "simulation_data_16_4.0.data"
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
        if sys.argv[i] == "-save":
            plot_file = sys.argv[i+1]
    #
    # plot puncture data...
    #
    puncture_plot(data_file, plot_file)


def usage():
    print("Creates plots from data files produced with simulation.py.")
    print("-"*30)
    print("The following options can be used to over-write default parameters")
    print(
        "\t-data: provide name of data_file [Default: simulation_data_16_4.0.data] ")
    print(
        "\t-save: if provided, plot will be saved to this file [Default: None]")
    print("-"*30)
    print("For example, to make a plot of the file simulation_data_24_6.0.data")
    print("and save to file example.pdf, call")
    print("\tpython plot_simulation.py -data simulation_data_24_6.0.data -save example.pdf")


if __name__ == '__main__':
    main()
