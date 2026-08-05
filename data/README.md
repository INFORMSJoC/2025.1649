# Benchmark Instances

This directory contains all benchmark instances used in the computational experiments reported in the paper.

The same set of benchmark instances is also provided in `src/circle_covering/data`, where they are directly accessed by the C++ implementation.

## File Format

Each instance is stored as a plain text file describing a polygonal region and, optionally, one or more holes.

The first section defines the outer boundary of the region:

```
C m
x1 y1
x2 y2
...
xm ym
```

where

* `C` indicates the container (outer polygon);
* `m` is the number of vertices of the container;
* `(xi, yi)` are the coordinates of the polygon vertices listed in order (clockwise or counterclockwise).

If the region contains holes, each hole is described by a separate section:

```
H k
x1 y1
x2 y2
...
xk yk
```

where

* `H` indicates a hole;
* `k` is the number of vertices of the hole;
* `(xi, yi)` are the coordinates of the hole boundary listed in order.

An instance may contain any number of hole sections. If no `H` section appears, the region contains no holes.

## Example

```
C 9
0.5 0.5
3.5 3
4 0.5
6 2
11 0.5
11 8
7 3
7.5 7
2 5
H 7
5 2.5
7 2.4
...
H 17
8.5 2
9 1.7
...
```

This example represents

* an outer polygon with **9** vertices;
* the first hole with **7** vertices;
* the second hole with **17** vertices.

## Notes

* Polygon vertices should be listed sequentially along each boundary.
* Multiple `H` sections indicate multiple holes within the same region.
* All benchmark instances in this directory are identical to those in `src/Circle Covering/data`, which are used directly by the implementation.
