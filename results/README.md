# Results

This directory contains the best solutions obtained in the computational experiments reported in the paper.

## Directory Structure

The results are organized according to the problem type and benchmark instance.

```
results/
├── constrained CRICMR/
│   ├── C11/
│   ├── C14H2/
│   ├── ...
│
└── unconstrained CRICMR/
    ├── Square/
    ├── Triangle/
    ├── ...
```

Each benchmark instance contains one or more solution files corresponding to different numbers of circles.

For example,

```
results/unconstrained CRICMR/Square/50.txt
```

stores the best solution found for covering the **Square** instance with **50 identical circles**.

## File Format

Each solution file is a plain text file with the following format:

```
n  r
x1  y1
x2  y2
...
xn  yn
```

where

* `n` is the number of circles;
* `r` is the optimized common covering radius;
* `(xi, yi)` denotes the center coordinates of the *i*-th circle.

## Solution Verification

A solution can be considered valid if it satisfies both the coverage requirement and the geometric constraints (for constrained instances).

When evaluating a solution, the following quantities should be checked:

* **Coverage violation (sum of squared errors):** For each Voronoi vertex of a circle's cell, if its distance to the corresponding circle center exceeds the covering radius `r`, the squared difference is accumulated. A value close to zero indicates that the region is completely covered.
* **Constraint violation (sum of squared errors):** (for constrained instances only) If a circle center lies outside the container or inside a hole, the squared distance to the nearest feasible boundary is accumulated. A value close to zero indicates that all circle centers satisfy the geometric constraints.

A solution is regarded as feasible when both quantities are close to zero (up to numerical precision).


## Notes

The files in this directory correspond to the best solutions reported in the paper and can be used for reproducing the computational results or for visualization and further analysis.
