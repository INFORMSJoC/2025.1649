# Benchmark Instances

This directory contains all the regions used to generate the benchmark instances which were used in our computational experiments.

The same set of benchmark regions is also provided in `src/circle_covering/data`.

## File Format

Each benchmark region is stored as a plain text file describing a polygonal region possibly with one or more holes.

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
* `(xi, yi)` are the coordinates of the vertices of polygonal container, which are listed in a counterclockwise order.

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
* `(xi, yi)` are the coordinates of the vertices of polygonal hole , which are listed in a counterclockwise order.

A polygonal region may contain holes. If no `H` appears, the region contains no holes.

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

* The vertices of each polygon should be listed in counterclockwise order.
* The region contains multiple holes if 'H' appears multiple times.
* All benchmark regions in this directory are the same as those in `src/circle_covering/data`.
