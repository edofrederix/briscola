# Chapter 8: Testing and validation

[Back to the table of contents](./0_start.md)
or [Previous chapter: Linear system solvers](./7_solvers.md)

## Unit tests

Briscola contains a number of semi-automated unit tests. These can be found in
`applications/test`. The test applications contained in each subfolder are
generally meant to test a number of functions related to a specific class. For
example, we can take a look at `Test-metrics.C` in `applications/test/metrics`.
In this test, a mesh is constructed from a `briscolaMeshDict`. The corresponding
mesh metrics are then systematically checked for both the colocated and the
staggered mesh. If any of the metrics are incorrect, an error will be given.

To run the test, let us change to the `applications/test/metrics` directory.
The test application can be compiled by running the build script:
```
./build.sh
```
To run the test, the run script can be used:
```
./run.sh
```
The test application will run, and return the following message if the test did
not encounter any errors: `Test metrics succesful`. If an error is encountered,
the specific test number where the error was encountered will be included in the
error message. With this test number, it is easy to find within `Test-metrics.C`
where the test went wrong.

Many other tests are also included with Briscola. Take a look at the `tests`
directory to see them. All of Briscola's test cases can be compiled and run in
the same way as described above. Alternatively, a script is provided to
automatically compile and run all tests at once. To do this, change to the
`applications/test` directory and use the run script:
```
./run.sh
```
This script will by default run up to 16 tests in parallel until all tests have
been run. An argument can also be given to the script to run a different number
of tests in parallel. All test directories can be cleaned and returned to their
original state with the `clean.sh` script.

Unit tests are an important quality control tool used throughout Briscola's
development, as they allow for the early identification of bugs when changes
are made to the code. Therefore, for users who intend to make changes or
additions to the code, it is recommended to run the test cases after modifying
the code to ensure that nothing was accidentally broken.

## Validation cases

Besides unit tests, a number of validation cases are also provided for Briscola.
These validation cases can serve as integral tests for different flow solvers.
In addition, they can highlight the impact of different solver choices (e.g.,
staggered vs. colocated) and models (e.g., different interface normal schemes in
two-phase simulations). As an example, let us have a look at the Hysing
validation case in `validation/Hysing`. This case is based on the same benchmark
which is also used in the [tutorial](./1_tutorial.md), so have a look there
first for more details on the case.

### Running the case

In the Hysing directory, a `template` folder contains a generic case folder for
the Hysing case, where a number of input files are parametrized such that
different inputs can be given using `m4`. These inputs are set by the
`Hysing.sh` script. If we take a look at this script, we can see that the four
different two-phase solvers are listed to be tested. Also different mesh
resolutions, parallel decompositions, pressure solvers, and interface normal and
curvature schemes are listed. Further down, we can see that the script loops
through all possible combinations of all of the parameters:
```
for I in "${!SOLVERS[@]}"; do
for J in "${!MESHES[@]}"; do
for K in "${!NPROCSPERBRICKSIDE[@]}"; do
for L in "${!PSOLVERS[@]}"; do
for M in "${!NORMALSCHEMES[@]}"; do
for N in "${!CURVATURESCHEMES[@]}"; do
    ...
```
then runs the case with the specified solver:
```
if [ "$NPROC" == "1" ]; then

    $SOLVER > log.$SOLVER

else

    mpirun \
        --bind-to none \
        --oversubscribe \
        -n $NPROC \
        $SOLVER -parallel > log.$SOLVER

fi
```
and finally runs a post-processing python script which creates plots to compare
the simulation results with the benchmark reference, as well as computes some
mean errors:
```
$PYTHON post.py log.$SOLVER
```
Since the `Hysing.sh` script launches a large number of smulations, it is
recommended to run it on a sufficient number of cores. At least 16 cores are
needed for some of the decompositions, but the script can also be run on more
cores using slurm:
```
sbatch Hysing.job
```
This will allocate 64 cores to the script, and a new simulation will be started
automatically if enough cores are free, untill all simulations have completed.
The number of cores can be also be increased or decreased in `Hysing.job`.

### Viewing the results

After all of the validation simulations are finished running, a few different
outputs can be inspected to see the outcomes. First, we can take a look at all
of the individual validation runs in the `runs` directory. The name of each
subdirectory here indicates the combination of two-phase solver, mesh resolution,
number of processors, pressure solver, normal scheme and curvature scheme.
Within each subdirectory, the evolution of the velocity and position of the
bubble is plotted and compared with the reference.

Looking at these plots for each of the individual validation runs can be tedious.
Therefore, the results are also summarized in `results.csv` within the main
Hysing directory. For each run, this file shows the mean error in bubble position
and in bubble velocity (`error 1 [%]` and `error 2 [%]` respectively). The
entries for `test 1` and `test 2` indicate whether these errors are within an
acceptable range. Additionally, the number of time steps and number of pressure
iterations are given for each run as well.

### Other cases

At present, four validation cases are available in the `validation` folder: the
lid-driven cavity (both colocated and staggered), the case of flow over a
cylinder with periodic vortex shedding, the Hysing case, and the laminar
Poiseuille flow case (both colocated and staggered, and with standard boundary
conditions as well as immersed boundary conditions). All validation cases can be
run at once using the `submitJobs.sh` script. It should be noted that running
this script will launch a large number of simulations, which may take a long
time and a lot of computational resources to run. All validation directories
can be cleaned and returned to their original state with the `clean.sh` script.

[Back to the table of contents](./0_start.md)
or [Next chapter: Current limitations](./9_limitations.md)

