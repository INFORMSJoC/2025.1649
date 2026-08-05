[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# A novel and efficient formulation and solution approach for the region coverage problem with identical circles

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the MIT License.

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on the paper
[A novel and efficient formulation and solution approach for the region coverage problem with identical circles](https://doi.org/10.1287/ijoc.2025.1649) by Q.H. Wu, M.X. Li, X.J. Lai, J.K. Hao, and Y.Z. Qiu.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2025.1649

https://doi.org/10.1287/ijoc.2025.1649.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{CircleCovering2026,
  author =        {Qinghua Wu, Minxuan Li, Xiangjing Lai, Jin-Kao Hao and Yuzhuo Qiu},
  publisher =     {INFORMS Journal on Computing},
  title =         {A novel and efficient formulation and solution approach for the region coverage problem with identical circles},
  year =          {2026},
  doi =           {10.1287/ijoc.2025.1649.cd},
  url =           {https://github.com/INFORMSJoC/2025.1649},
  note =          {Available for download at https://github.com/INFORMSJoC/2025.1649},
}
```

## Description

This repository contains a C++ implementation of the algorithm proposed in the paper for solving the **region coverage problem with identical circles**. Given a polygonal region (possibly containing holes) and a specified number of circles, the program searches for circle placements that minimize the common covering radius.

The repository also includes the benchmark instances, experimental results, and auxiliary scripts used in the computational study reported in the paper.

## Running the Program

To generate the executable code, compile the project using CMake in the `src/circle_covering` directory.

### Usage

```bash
./CircleCovering Region n Constraint
```

* `Region` denotes the benchmark instance (e.g., `Square`, `Triangle`, `C14H2`, `C101H3`).
* `n` is the number of circles (positive integer).
* `Constraint` specifies the problem type:

  * `0`: unconstrained region coverage problem;
  * `1`: constrained region coverage problem.

Example:

```bash
./CircleCovering C14H2 50 1
```

If no command-line arguments are provided, the program will request the above information interactively.

*Note:* Example batch scripts for Linux clusters are provided in the [scripts](scripts) directory.

## Materials

This repository includes the following materials:

* *Benchmark instances for both constrained and unconstrained region coverage problems.* (See the [data](data) directory for the details.)
* *C++ source code of the proposed circle covering algorithm.* (See [the source codes](src/circle_covering) directory for the details.)
* *Example scripts for compiling and running the program on Linux clusters.* (See the [scripts](scripts) directory for the details.)
* *Best solutions obtained in the computational experiments reported in the paper.* (See the [results](results) directory for the details.)
* *Supplementary materials accompanying the paper, including additional experimental data and documentation.* (See the [docs](docs) directory for the details.)

**Note:** Each subdirectory contains its own README file (when applicable) describing the corresponding file formats and contents.