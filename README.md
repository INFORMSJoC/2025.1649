[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# A novel and efficient formulation and solution approach for the region coverage problem with identical circles

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

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

This repository provides a C++ implementation of the algorithm proposed in the paper **A novel and efficient formulation and solution approach for the region coverage problem with identical circles** by Qinghua Wu, Minxuan Li, Xiangjing Lai, Jin-Kao Hao, and Yuzhuo Qiu.

The repository includes the source code of algorithm, the benchmark instances,  the best solutions found in the study,  the second part of online supplement, and the script to compile the source code of algorithm.  

## Running the Program

To generate the executable code of algorithm, compile the project in the`src/circle_covering` directory using CMake.

### Usage

```bash
./CircleCovering Region n Constraint
```

* `Region` denotes the benchmark region (e.g., `Square`, `Triangle`, `C14H2`, `C101H3`).
* `n` is a positive integer indicating the number of circles.
* `Constraint` is a binary integer flag  used to select the problem type:

  * `0` for the unconstrained region coverage problem;
  * `1` for the constrained region coverage problem.

Example:

```bash
./CircleCovering C14H2 50 1
```

## Materials

This repository includes the following materials:

* *Benchmark regions for both constrained and unconstrained region coverage problems.* (See the [data](data) directory for the details.)
* *C++ source code of the proposed  algorithm.* (See the [src](src/circle_covering) directory for the details.)
* *Script for compiling and running the program.* (See the [scripts](scripts) directory for the details.)
* *Best solutions found in the paper.* (See the [results](results) directory for the details.)
* *The second part of online supplement of the paper.* (See the [docs](docs) directory for the details.)

**Note:** Each subdirectory contains its own README file.
