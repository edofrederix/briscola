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
TASKFILE="/tmp/tasks.$$"

# Simulation parameters

MESHES=(uniform graded edgeGraded)
NBRICKS=(1 4 4)
NPROCSPERBRICKSIDE=(1 2)
LSCHEMES=(linearGauss)
DSCHEMES=(linearGauss midPointGauss)
RES=(400 1000)
COARSEMODES=(smooth direct)
RKSCHEMES=(forwardEuler RK3 Ascher222 CNAB)

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
    "laplacian scheme," \
    "divergence scheme," \
    "Re," \
    "coarse mode," \
    "Time scheme," \
    "error 1 [%]," \
    "error 2 [%]," \
    "test 1," \
    "test 2," \
    "number of time steps," \
    "mean number of pressure iters" > $CURR/$CSV

rm -f $TASKFILE.[0-9]+

getNumTasks () {

    COUNT=0

    for PID in $(jobs -p -r); do

        while [ ! -f "$TASKFILE.$PID" ]; do

            sleep 1

        done

        COUNT=$(($COUNT + $(cat $TASKFILE.$PID)))

    done

    echo $COUNT
}

# Handle SIGTERM

PIDS=()

cleanup() {
    echo "Caught signal, killing children..."
    kill -TERM "${PIDS[@]}" 2>/dev/null
    wait
    exit 1
}

trap cleanup SIGINT SIGTERM

##

for I in "${!MESHES[@]}"; do
for J in "${!NPROCSPERBRICKSIDE[@]}"; do
for K in "${!LSCHEMES[@]}"; do
for L in "${!DSCHEMES[@]}"; do
for M in "${!RES[@]}"; do
for N in "${!COARSEMODES[@]}"; do
for O in "${!RKSCHEMES[@]}"; do

    sleep 1

    MESH=${MESHES[$I]}
    NBRICK=${NBRICKS[$I]}
    NPROCPERBRICKSIDE=${NPROCSPERBRICKSIDE[$J]}
    LSCHEME=${LSCHEMES[$K]}
    DSCHEME=${DSCHEMES[$L]}
    RE=${RES[$M]}
    COARSEMODE=${COARSEMODES[$N]}
    RKSCHEME=${RKSCHEMES[$O]}

    NPROCX=$NPROCPERBRICKSIDE
    NPROCY=$NPROCPERBRICKSIDE

    NPROC=$(echo "print(int($NPROCX*$NPROCY*$NBRICK))" | python)

    NU=$(echo "print(1.0/$RE)" | python)

    CASE="$MESH-$NPROC-$LSCHEME-$DSCHEME-$RE-$COARSEMODE-$RKSCHEME"

    while [ "$(($(getNumTasks) + $NPROC))" -gt $NTASKS ]; do

        echo \
            Procs running = $(getNumTasks), \
            Procs needed = $NPROC

        sleep 1

    done

    (
        echo "Starting $CASE"

        echo $NPROC > $TASKFILE.$BASHPID

        cp -r $TEMPLATE $RUNDIR/$CASE

        cd $RUNDIR/$CASE

        ./prep.sh $MESH $NPROCX $NPROCY $LSCHEME $DSCHEME $NU $COARSEMODE $RKSCHEME

        if [ "$RKSCHEME" == "CNAB" ]; then

            SOLVER2=${SOLVER}CNAB

        else

            SOLVER2=$SOLVER

        fi

        if [ "$NPROC" == "1" ]; then

            $SOLVER2 > log.$SOLVER

        else

            mpirun \
                --bind-to none \
                --oversubscribe \
                -n $NPROC \
                $SOLVER2 -parallel > log.$SOLVER

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

    ) & PIDS+=($!)

done
done
done
done
done
done
done

wait

rm -f $TASKFILE.[0-9]+
