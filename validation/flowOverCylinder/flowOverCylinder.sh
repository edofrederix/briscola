#!/bin/bash

# By default 16 cores are used. Simulations are sequentially started only when
# resources are available. Output is written to results.csv.

NTASKS="${SLURM_NTASKS:-16}"

# Some parameters

RUNDIR="runs"
CSV="results.csv"
TEMPLATE="template"
SOLVER="briscolaColocated"
PYTHON="python3"
NPROC=12

# Simulation parameters

MESHES=(6 8 12)
DSCHEMES=(linearGauss midPointGauss)
GRADSCHEMES=(linearGauss midPointGauss)
RES=(200 400)
RKSCHEMES=(backwardEuler RK3 Ascher222 CNAB)

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

##

cd $TEMPLATE
wmake -s code
cd $CURR

##


echo \
    "Mesh," \
    "divergence scheme," \
    "gradient scheme," \
    "Re," \
    "Time scheme," \
    "error [%]," \
    "test," \
    "number of time steps," \
    "mean number of pressure iters" > $CURR/$CSV

##

source ../procManagement.sh

for I in "${!MESHES[@]}"; do
for J in "${!DSCHEMES[@]}"; do
for K in "${!GRADSCHEMES[@]}"; do
for L in "${!RES[@]}"; do
for M in "${!RKSCHEMES[@]}"; do

    sleep 1

    MESH=${MESHES[$I]}
    DSCHEME=${DSCHEMES[$J]}
    GRADSCHEME=${GRADSCHEMES[$K]}
    RE=${RES[$L]}
    RKSCHEME=${RKSCHEMES[$M]}

    NU=$(echo "print(1.0/$RE)" | python)

    CASE="$MESH-$DSCHEME-$GRADSCHEME-$RE-$RKSCHEME"

    wait_for_procs $NPROC $NTASKS

    (
        echo "Starting $CASE"

        cp -r $TEMPLATE $RUNDIR/$CASE

        cd $RUNDIR/$CASE

        ./prep.sh $MESH $DSCHEME $GRADSCHEME $NU $RKSCHEME

        if [ "$RKSCHEME" == "CNAB" ]; then

            SOLVER2=${SOLVER}CNAB

        else

            SOLVER2=$SOLVER

        fi

        if [ "$NPROC" == "1" ]; then
            srun --exclusive -n $NPROC $SOLVER2 > log.$SOLVER
        else
            srun --exclusive -n $NPROC $SOLVER2 -parallel > log.$SOLVER
        fi

        ##

        $PYTHON post.py log.$SOLVER > result.txt

        E=$(sed -n '1p' < result.txt)
        P=$(sed -n '2p' < result.txt)

        NDT=$(sed -n '3p' < result.txt)
        NITER=$(sed -n '4p' < result.txt)

        echo \
            "$MESH," \
            "$DSCHEME," \
            "$GRADSCHEME," \
            "$RE," \
            "$RKSCHEME," \
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

wait
