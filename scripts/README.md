# Cluster / Batch Submission (HPC)
Note: the CMake target/executable is named `CircleCovering` as defined in `CMakeLists.txt` (see `add_executable`); use that name when running or editing `run.sh`.

When testing on a Linux cluster using batch job submission (for example Slurm), you may prefer non-interactive runs. For that case, the repository contains an alternative `main(int argc, char *argv[])` (kept in the source as a commented example) which accepts `Region n CONSTRAINT` as command-line arguments.

You can use or adapt the provided `run.sh` batch script as an example. `run.sh` builds the project and launches multiple non-interactive runs (backgrounded). 

Minimal example (inside a Slurm script or bash job):

```bash
cd build
cmake .. && make -j
./CircleCovering C4 100 1 &
./CircleCovering C92H4 50 1 &
wait
```

The provided `run.sh` in the repository is a ready-made Slurm batch script that demonstrates module loading, building, and multiple runs. Edit it to match your cluster environment and executable name.