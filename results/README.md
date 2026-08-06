# Results

This directory contains the best solutions obtained in our computational experiments.

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

Each innermost subfolder contains the solution files of several benchmark instances sharing the same region.

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


## Notes

Each innermost subfolder of this directory corresponds to a benchmark region, and contains the coordinate file of benchmark region and solution files of several instances with different numbers of circles.
