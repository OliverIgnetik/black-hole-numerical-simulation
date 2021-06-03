# Numerical Simulation Black Hole

## Overview

This is an attempt at modelling initial data and the associated extrinsic curvature of a black hole with linear momentum. The associated extrinsic curvature is calculated in close proximity to the equatorial plane of the black hole. The approach is based on the puncture method outlined in Baumgarte et al. [3].

![Extrinsic Curvature of Spacetime manifold](img/bh_nasa.gif)

Image credit : [NASA](https://www.nasa.gov/feature/goddard/2019/nasa-visualization-shows-a-black-hole-s-warped-world/)

## Motivation

This project was an attempt to make use of the understanding gained in recent self-directed study of general relativity, the associated tensor calculus and differential geometry.

## Usage

The script `simulation.py` makes use of the puncture method as described in to construct initial data describing black holes.
`simulation.py` has a number of flags that can be used to customize the output of the script:

### `simulation.py flags`

```
-grid_dim (grid points in each dimensions)

-boundary (location of outer boundary)

-loc_x (x coordinate of black hole)

-loc_y (y coordinate of black hole)

-loc_z (z coordinate of black hole)

-p_x (x momentum of black hole)

-p_y (y momentum of black hole)

-p_z (z momentum of black hole)

-tol (tolerance for elliptic solver)

-it_max (maximum number of iterations)
```

The environment used for this code can be easily set up using the `anaconda-project.yml` file.
The script will will produce a file called `simulation_data_{grid_dim}_{x_boundary}.data`.

The data file contains values of the function u dependent on x and y on a plane of constant z in close proximity to the equatorial plane.

Once the data has been constructed using `simulation.py` we can make use of `plot_simulation.py` for visualization.

### `plot_simulation.py flags`

```
-data (puncture data filename)
-save (save filename for image)
```

Below is an example of initial data of a black hole with dimensionless linear momentum P = (1,0,0) located at the origin.

![Extrinsic Curvature of Spacetime manifold](img/example_plot.png)

## Packages

Please see `anaconda-project.yml` for environment dependencies and package versions

- `python 3.95`
- `numpy 1.20.2`
- `scipy 1.6.2`
- `matplotlib 3.3.4`

## Features

- Constructs initial data for a black hole with linear momentum using the puncture method
- Configurable boundary and number of grid points used in simulation
- Configurable linear momentum and location of the black hole
- `elliptic_solver.py` interface with `scipy.la` for solving constraint equations

## Future Work

- Write up an accompanying latex formatted article explaining the mathematics of the simulation
- Generalize this approach to the moving puncture method to model binary neutron star precession
- Add a Graphical User Interface for easier manipulation of parameters

## Non-Exhaustive References

### Numerical Relativity

[1]E. Gourgoulhon, _3+1 formalism in general relativity_. Heidelberg: Springer, 2012.

[2]M. Shibata, _Numerical relativity_. Singapore: World Scientific Publishing Co. Pte. Ltd., 2016.

[3]T. Baumgarte and S. Shapiro, _Numerical relativity_. Cambridge: Cambridge University Press, 2010.

### General Relativity and Differential Geometry

[4] B. Schutz, _A first course in general relativity_. Cambridge [etc.]: Cambridge University Press, 2018.

[5] C. Misner, K. Thorne and J. Wheeler, _Gravitation_. New York: Freeman, 1995.

[6] S. Chandrasekhar, _The mathematical theory of black holes_. Oxford u.a.: Clarendon Press, 1992.

[7] D. Neuenschwander, _Tensor calculus for physics_. John Hopkins University Press, 2015.

[8] S. Carroll, _Spacetime and geometry_. Harlow: Pearson, 2014.

[9] N. Straumann, _General relativity_. Heidelberg: Springer, 2013.

[10] J. Hubbard and B. Hubbard, _Vector calculus, linear algebra, and differential forms_, 5th ed. New York: Matrix Editions, 2015.
