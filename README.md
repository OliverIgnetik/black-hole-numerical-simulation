# Numerical Simulation Black Hole

## Summary

This is an attempt at modelling initial data and the associated conformal deviation of spacetime around a black hole with linear momentum. The curvature is calculated in close proximity to the equatorial plane of the black hole. The approach is based on the puncture method outlined in _Shibata_ [2].

![Extrinsic Curvature of Spacetime manifold](img/bh_nasa.gif)

Image credit : [NASA](https://www.nasa.gov/feature/goddard/2019/nasa-visualization-shows-a-black-hole-s-warped-world/)

<!-- Please see [numerical_simulation_theory.pdf]() included in this repository for a detailed explanation of:

- Decompostion of Einstein's Field Equations into constraint and evolution equations using 3+1 foliation
- Puncture method for solving the constraint equation -->

## Methodology Overview

This approach recasts Einstein's field equations using the well known 3+1 decomposition of spacetime. This decomposition splits the gravitational field equations into constraint and evolution equations. The constraint equations are crucial for imposing conditions on the gravitational fields at any moment in time.

By the use of a conformal decomposition of the the constraint equations we can describe the initial gravitational field of the black hole using the puncture method. Note that this method is generalizable to model the initial data of multiple black holes with both linear and angular momentum, this is because of the resulting linearity of the decoupled constraint equation. This means me can model find solutions for spacetimes more general then that described by the Schwarzschild metric; such as those in which black holes have angular and linear momentum. Despite the fact that this approach is valid for these more general solutions, we have not yet incoporated the angular momentum of the black hole in the simulation. The code presented in this repository solves the constraint equations for initial data using the puncture method for spherically symmetric black holes with _linear momentum_.

Please see either _Gourgoulhon_ [1] or _Shibata_ [2] for a more in depth explanation of both the mathematics and the 3+1 decomposition of spacetime.

<!-- Please see [numerical_simulation_theory.pdf]() for further clarity. -->

## Motivation

This project was an attempt to make use of the understanding gained in recent self-directed study of general relativity, tensor calculus and differential geometry.

## Usage

The script `simulation.py` solves the constraint equation by making use of the puncture method as described in Shibata [2], to construct initial data describing the conformal deviation of spacetime on the equatorial plane of a black hole. `simulation.py` has a number of flags that can be used to customize the output of the script:

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

### Example simulation output

Below is an example of initial data of a black hole with dimensionless linear momentum `P = (1,0,0)` located at the origin. `u` is the correction of the Schwarzschild conformal factor.

![Conformal deviation of Spacetime manifold](img/example_plot.png)

Note that `M` as shown in the graph above is the _puncture mass_ which dominates the conformal factor close to the black hole's center. We use the _puncture mass_ to make all variables such as coordinate position and momentum dimensionless.

## Packages

Please see `anaconda-project.yml` for environment dependencies and package versions

- `python 3.95`
- `numpy 1.20.2`
- `scipy 1.6.2`
- `matplotlib 3.3.4`

## Features

- Solves the constraint equation and constructs initial data for a black hole with linear momentum using the puncture method
- Configurable boundary and number of grid points used in simulation
- Configurable linear momentum and location of the black hole
- `elliptic_solver.py` interface with `scipy.la` for solving constraint equations

## Future Work

- Write up an accompanying latex formatted article explaining the mathematics of the simulation
- Generalize code to account for black holes with angular momentum
- Work on **incorporating the initial data in the evolution equations**
- Generalize this approach to the moving puncture method to model binary neutron star precession
- Add a Graphical User Interface for easier manipulation of parameters

## Non-Exhaustive References

### Numerical Relativity

[1] E. Gourgoulhon, _3+1 formalism in general relativity_. Heidelberg: Springer, 2012.

[2] M. Shibata, _Numerical relativity_. Singapore: World Scientific Publishing Co. Pte. Ltd., 2016.

[3] T. Baumgarte and S. Shapiro, _Numerical relativity_. Cambridge: Cambridge University Press, 2010.

### General Relativity and Differential Geometry

[4] B. Schutz, _A first course in general relativity_. Cambridge [etc.]: Cambridge University Press, 2018.

[5] C. Misner, K. Thorne and J. Wheeler, _Gravitation_. New York: Freeman, 1995.

[6] S. Chandrasekhar, _The mathematical theory of black holes_. Oxford u.a.: Clarendon Press, 1992.

[7] D. Neuenschwander, _Tensor calculus for physics_. John Hopkins University Press, 2015.

[8] S. Carroll, _Spacetime and geometry_. Harlow: Pearson, 2014.

[9] N. Straumann, _General relativity_. Heidelberg: Springer, 2013.

[10] J. Hubbard and B. Hubbard, _Vector calculus, linear algebra, and differential forms_, 5th ed. New York: Matrix Editions, 2015.
