#!/bin/bash

# By default 16 cores are used. Simulations are sequentially started only when
# resources are available. Output is written to results.csv.

NTASKS="${SLURM_NTASKS:-16}"

# Some parameters

RUNDIR="runs"
CSV="results.csv"
TEMPLATE="template"

# Simulation parameters

MESHES=(32 64)
MODES=(mapped penalization Vreman Fadlun Mittal)
NPROCSPERBRICKSIDE=(1 2 4)
RKSCHEMES=(RK3 Ascher222 CNAB)
SOLVERS=(briscolaColocated briscolaStaggered)
PSOLVERS=(MG FFT)

##

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

CURR=$(pwd)

if [ -d "$RUNDIR" ]; then
    rm -fr $RUNDIR
fi

mkdir -p $RUNDIR

echo \
    "Mesh," \
    "Nprocs," \
    "Mode," \
    "Pressure solver," \
    "Time scheme," \
    "Solver," \
    "Error [%]," \
    "Test," \
    "Number of time steps," \
    "Mean number of pressure iters" > $CURR/$CSV
##

source ../procManagement.sh

for I in "${!MESHES[@]}"; do
for J in "${!MODES[@]}"; do
for K in "${!NPROCSPERBRICKSIDE[@]}"; do
for L in "${!PSOLVERS[@]}"; do
for M in "${!RKSCHEMES[@]}"; do
for N in "${!SOLVERS[@]}"; do

    sleep 1

    MESH=${MESHES[$I]}
    MODE=${MODES[$J]}
    NPROCPERBRICKSIDE=${NPROCSPERBRICKSIDE[$K]}
    PSOLVER=${PSOLVERS[$L]}
    RKSCHEME=${RKSCHEMES[$M]}
    SOLVER=${SOLVERS[$N]}

    NPROCX=$NPROCPERBRICKSIDE
    NPROCY=$NPROCPERBRICKSIDE

    NPROC=$(echo "print(int($NPROCX*$NPROCY))" | python)

    CASE="$MESH-$NPROC-$MODE-$PSOLVER-$RKSCHEME-$SOLVER"

    wait_for_procs $NPROC $NTASKS

    (
        echo "Starting $CASE"

        cp -r $TEMPLATE $RUNDIR/$CASE

        cd $RUNDIR/$CASE

        ./prep.sh $MESH $NPROCX $NPROCY $MODE $PSOLVER $RKSCHEME

        if [ "$NPROC" == "1" ]; then
            srun --exclusive -n $NPROC -n 1 $SOLVER > log.$SOLVER
        else
            srun --exclusive -n $NPROC $SOLVER -parallel > log.$SOLVER
        fi

        ##

        $PYTHON post.py log.$SOLVER

        E=0
        P=failed
        NDT=0
        NITER=0

        if [ -f "result.txt" ]; then

            E=$(sed -n '1p' < result.txt)
            P=$(sed -n '2p' < result.txt)

            NDT=$(sed -n '3p' < result.txt)
            NITER=$(sed -n '4p' < result.txt)

        fi

        echo \
            "$MESH," \
            "$NPROC," \
            "$MODE," \
            "$PSOLVER," \
            "$RKSCHEME," \
            "$SOLVER," \
            "$E," \
            "$P," \
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

wait
