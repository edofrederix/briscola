#!/bin/bash

# By default 16 cores are used. Simulations are sequentially started only when
# resources are available. Output is written to results.csv.

NTASKS="${SLURM_NTASKS:-16}"

# Some parameters

RUNDIR="runs"
CSV="results.csv"
TEMPLATE="template"
PYTHON="python3"

# Simulation parameters

SOLVERS=(briscolaColocatedTwoPhase briscolaStaggeredTwoPhase)
MESHES=(32 64)
NPROCSPERBRICKSIDE=(2 4)
PSOLVERS=(MG split)
NORMALSCHEMES=(MYC LSGIR)
CURVATURESCHEMES=(SHF CV)
RKSCHEMES=(Ascher222 CNABM DIRK2)

##

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

##

CURR=$(pwd)

##

if [ -d "$RUNDIR" ]; then
    rm -fr $RUNDIR
fi

mkdir -p $RUNDIR

echo \
    "solver," \
    "mesh," \
    "Nprocs," \
    "pressure solver," \
    "normal scheme," \
    "curvature scheme," \
    "Runge-Kutta scheme," \
    "error 1 [%]," \
    "error 2 [%]," \
    "test 1," \
    "test 2," \
    "number of time steps," \
    "mean number of pressure iters" > $CURR/$CSV

##

source ../procManagement.sh

for I in "${!SOLVERS[@]}"; do
for J in "${!MESHES[@]}"; do
for K in "${!NPROCSPERBRICKSIDE[@]}"; do
for L in "${!PSOLVERS[@]}"; do
for M in "${!NORMALSCHEMES[@]}"; do
for N in "${!CURVATURESCHEMES[@]}"; do
for O in "${!RKSCHEMES[@]}"; do

    sleep 1

    SOLVER=${SOLVERS[$I]}
    MESH=${MESHES[$J]}
    NPROCPERBRICKSIDE=${NPROCSPERBRICKSIDE[$K]}
    PSOLVER=${PSOLVERS[$L]}
    NORMALSCHEME=${NORMALSCHEMES[$M]}
    CURVATURESCHEME=${CURVATURESCHEMES[$N]}
    RKSCHEME=${RKSCHEMES[$O]}

    NPROCX=$NPROCPERBRICKSIDE
    NPROCY=$NPROCPERBRICKSIDE

    NPROC=$(echo "print(int($NPROCX*$NPROCY))" | python)

    MESHX=$MESH
    MESHY=$(echo "print(int(2*$MESH))" | python)

    CASE="$SOLVER-$MESH-$NPROC-$PSOLVER-$NORMALSCHEME-$CURVATURESCHEME-$RKSCHEME"

    wait_for_procs $NPROC $NTASKS

    (
        echo "Starting $CASE"

        cp -r $TEMPLATE $RUNDIR/$CASE

        cd $RUNDIR/$CASE

        ./prep.sh \
            $SOLVER \
            $MESHX \
            $MESHY \
            $NPROCX \
            $NPROCY \
            $PSOLVER \
            $NORMALSCHEME \
            $CURVATURESCHEME \
            $RKSCHEME

        if [ "$NPROC" == "1" ]; then
            srun --exclusive -n $NPROC $SOLVER > log.$SOLVER
        else
            srun --exclusive -n $NPROC $SOLVER -parallel > log.$SOLVER
        fi

        ##

        $PYTHON post.py log.$SOLVER

        if [ -f "result.txt" ]; then

            E1=$(sed -n '1p' < result.txt)
            E2=$(sed -n '2p' < result.txt)

            P1=$(sed -n '3p' < result.txt)
            P2=$(sed -n '4p' < result.txt)

            NDT=$(sed -n '5p' < result.txt)
            NITER=$(sed -n '6p' < result.txt)

        fi

        echo \
            "$SOLVER," \
            "$MESH," \
            "$NPROC," \
            "$PSOLVER," \
            "$NORMALSCHEME," \
            "$CURVATURESCHEME," \
            "$RKSCHEME," \
            "$E1," \
            "$E2," \
            "$P1," \
            "$P2," \
            "$NDT," \
            "$NITER" \ >> $CURR/$CSV

        ##

        cd $CURR

        echo "Finished $CASE"

    ) &

    store_procs $NPROC $!

done
done
done
done
done
done
done

wait
