# acoustic_knots: an exploration of how knots can arise in acoustic signatures

Copyright (c) 2020, Michael Robinson

## Synopsis

If an object is illuminated by sound as it is rotated, the echos will recur once the object has been turned through 360 degrees.  If the object is symmetric, then the echos will recur sooner.  If the object can be decomposed into several symmetric parts, as suggested by a Fourier series, then the echos won't necessarily recur as a whole until the object is rotated 360 degrees, but the decomposition can still be discerned.

Thinking of each position as a point, we might take the echos as specifying coordinates for that point.  The trajectory through "echo space" is somehow characteristic of the object.  The Fourier decomposition suggests that the trajectory can be thought of as a torus knot.  Topological classification of torus knots works if the torus is 2-dimensional, but not if it is higher dimensional.  (If the torus is 3-dimensional, then all torus knots are isotopic to the unknot!)

This repository contains a small acoustic simulator and some visualizations of the data as torus knots.

## Manifest

* `isotropicMatrix`: Isotropic Green's function for the linear scalar wave equation.  Translation: propagate a wave from a source point to a sample point, where it is measured.   

* `knotify`: Compute the residual from assuming a signal is a 2-dimensional torus knot of a given knot type.

* `toroidal_process`: Main script for the experiment.

* `pcoa`: Compute principal components for a matrix of data.  (Standard, but I didn't want any dependencies!)
