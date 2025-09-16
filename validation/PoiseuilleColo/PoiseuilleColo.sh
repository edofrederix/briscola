#!/bin/bash

# By default 16 cores are used. Simulations are sequentially started only when
# resources are available. Output is written to results.csv.

NTASKS="${SLURM_NTASKS:-16}"

# Some parameters

RUNDIR="runs"
CSV="results.csv"
TEMPLATE="template"
SOLVER="briscolaColocated"
TASKFILE="/tmp/tasks.$$"

# Simulation parameters

MESHES=(32 64)
IBMS=(none penalization Vreman Fadlun Mittal)
NPROCSPERBRICKSIDE=(1 2 4)
PSOLVERS=(MG FFT Krylov)
RKSCHEMES=(RK3 Ascher222 CNAB)

##

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

CURR=$(pwd)

cd $TEMPLATE
wmake -silent code
cd $CURR

if [ -d "$RUNDIR" ]; then
    rm -fr $RUNDIR
fi

mkdir -p $RUNDIR

echo \
    "Mesh," \
    "Nprocs," \
    "IBM," \
    "pressure solver," \
    "Time scheme," \
    "error [%]," \
    "test," \
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
for J in "${!IBMS[@]}"; do
for K in "${!NPROCSPERBRICKSIDE[@]}"; do
for L in "${!PSOLVERS[@]}"; do
for M in "${!RKSCHEMES[@]}"; do

    sleep 1

    MESH=${MESHES[$I]}
    IBM=${IBMS[$J]}
    NPROCPERBRICKSIDE=${NPROCSPERBRICKSIDE[$K]}
    PSOLVER=${PSOLVERS[$L]}
    RKSCHEME=${RKSCHEMES[$M]}

    NPROCX=$NPROCPERBRICKSIDE
    NPROCY=$NPROCPERBRICKSIDE

    NPROC=$(echo "print(int($NPROCX*$NPROCY))" | python)

    CASE="$MESH-$NPROC-$IBM-$PSOLVER-$RKSCHEME"

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

        ./prep.sh $MESH $NPROCX $NPROCY $IBM $PSOLVER $RKSCHEME

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

        $PYTHON post.py log.$SOLVER

        E=$(sed -n '1p' < result.txt)
        P=$(sed -n '2p' < result.txt)

        NDT=$(sed -n '3p' < result.txt)
        NITER=$(sed -n '4p' < result.txt)

        echo \
            "$MESH," \
            "$NPROC," \
            "$IBM," \
            "$PSOLVER," \
            "$RKSCHEME," \
            "$E," \
            "$P," \
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
