#!/bin/bash

# By default 16 cores are used. Simulations are sequentially started only when
# resources are available. Output is written to results.csv.

NTASKS="${SLURM_NTASKS:-16}"

# Some parameters

RUNDIR="runs"
CSV="results.csv"
TEMPLATE="template"
SOLVER="briscolaStaggered"
PYTHON="python3"

# Simulation parameters

MESHES=(uniform graded)
NBRICKS=(1 4)
NPROCSPERBRICKSIDE=(1 2)
LSCHEMES=(linearGauss)
DSCHEMES=(linearGauss midPointGauss)
RES=(100 400 1000)
COARSEMODES=(smooth direct)

##

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

CURR=$(pwd)

##

if [ -d "$RUNDIR" ]; then
    rm -fr $RUNDIR
fi

mkdir -p $RUNDIR

echo \
    "Mesh," \
    "Nprocs," \
    "laplacian scheme," \
    "divergence scheme," \
    "Re," \
    "coarse mode," \
    "error 1 [%]," \
    "error 2 [%]," \
    "test 1," \
    "test 2," \
    "number of time steps," \
    "mean number of pressure iters" > $CURR/$CSV

for I in "${!MESHES[@]}"; do
for J in "${!NPROCSPERBRICKSIDE[@]}"; do
for K in "${!LSCHEMES[@]}"; do
for L in "${!DSCHEMES[@]}"; do
for M in "${!RES[@]}"; do
for N in "${!COARSEMODES[@]}"; do

    sleep 1

    MESH=${MESHES[$I]}
    NBRICK=${NBRICKS[$I]}
    NPROCPERBRICKSIDE=${NPROCSPERBRICKSIDE[$J]}
    LSCHEME=${LSCHEMES[$K]}
    DSCHEME=${DSCHEMES[$L]}
    RE=${RES[$M]}
    COARSEMODE=${COARSEMODES[$N]}

    NPROCX=$NPROCPERBRICKSIDE
    NPROCY=$NPROCPERBRICKSIDE

    NPROC=$(echo "print(int($NPROCX*$NPROCY*$NBRICK))" | python)

    NU=$(echo "print(1.0/$RE)" | python)

    CASE="$MESH-$NPROC-$LSCHEME-$DSCHEME-$RE-$COARSEMODE"

    while [ "$(($(ps aux | grep $SOLVER | grep -v mpirun | grep -v grep | wc -l) + $NPROC))" \
        -gt $NTASKS ]; do

        sleep 1

    done

    (
        echo "Starting $CASE"

        cp -r $TEMPLATE $RUNDIR/$CASE

        cd $RUNDIR/$CASE

        ./prep.sh $MESH $NPROCX $NPROCY $LSCHEME $DSCHEME $NU $COARSEMODE

        if [ "$NPROC" == "1" ]; then

            $SOLVER > log.$SOLVER

        else

            mpirun \
                --bind-to none \
                --oversubscribe \
                -n $NPROC \
                $SOLVER -parallel > log.$SOLVER

        fi

        ##

        $PYTHON post.py $RE log.$SOLVER

        E1=$(sed -n '1p' < result.txt)
        E2=$(sed -n '2p' < result.txt)

        P1=$(sed -n '3p' < result.txt)
        P2=$(sed -n '4p' < result.txt)

        NDT=$(sed -n '5p' < result.txt)
        NITER=$(sed -n '6p' < result.txt)

        echo \
            "$MESH," \
            "$NPROC," \
            "$LSCHEME," \
            "$DSCHEME," \
            "$RE," \
            "$COARSEMODE," \
            "$E1," \
            "$E2," \
            "$P1," \
            "$P2," \
            "$NDT," \
            "$NITER" \ >> $CURR/$CSV

        ##

        echo "Finished $CASE"

    ) &

done
done
done
done
done
done

wait
