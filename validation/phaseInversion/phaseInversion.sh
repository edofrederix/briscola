#!/bin/bash

# By default 64 cores are used. Simulations are sequentially started only when
# resources are available. Output is written to results.csv.

NTASKS="${SLURM_NTASKS:-64}"

# Some parameters

RUNDIR="runs"
CSV="results.csv"
TEMPLATE="template"
PYTHON="python3"
TASKFILE="/tmp/tasks.$$"

# Simulation parameters

SOLVERS=(\
    briscolaColocatedTwoPhase \
    briscolaColocatedTwoPhaseCNAB \
    briscolaStaggeredTwoPhase \
    briscolaStaggeredTwoPhaseCNAB)
MESHES=(32 64 128)
NPROCSPERBRICKSIDE=(2 4)
PSOLVERS=(MG split)
NORMALSCHEMES=(Youngs MYC LSGIR)
CURVATURESCHEMES=(SHF CV)

##

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

##

CURR=$(pwd)

for I in "${!SOLVERS[@]}"; do

    SOLVER=${SOLVERS[$I]}

    cd $TEMPLATE
    wmake -silent code.$SOLVER

    cd $CURR

done

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
    "error 1 [%]," \
    "error 2 [%]," \
    "error 3 [%]," \
    "error 4 [%]," \
    "test 1," \
    "test 2," \
    "test 3," \
    "test 4," \
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

for I in "${!SOLVERS[@]}"; do
for J in "${!MESHES[@]}"; do
for K in "${!NPROCSPERBRICKSIDE[@]}"; do
for L in "${!PSOLVERS[@]}"; do
for M in "${!NORMALSCHEMES[@]}"; do
for N in "${!CURVATURESCHEMES[@]}"; do

    SOLVER=${SOLVERS[$I]}
    MESH=${MESHES[$J]}
    NPROCPERBRICKSIDE=${NPROCSPERBRICKSIDE[$K]}
    PSOLVER=${PSOLVERS[$L]}
    NORMALSCHEME=${NORMALSCHEMES[$M]}
    CURVATURESCHEME=${CURVATURESCHEMES[$N]}

    if ! { [ "$MESH" -eq 32 ] && [ "$NPROCPERBRICKSIDE" -eq 4 ]; } && \
       ! { [ "$MESH" -eq 128 ] && [ "$NPROCPERBRICKSIDE" -eq 2 ]; }; then

        sleep 1

        NPROCX=$NPROCPERBRICKSIDE

        NPROC=$(echo "print(int($NPROCX**3))" | python)

        MESHX=$MESH

        CASE="$SOLVER-$MESH-$NPROC-$PSOLVER-$NORMALSCHEME-$CURVATURESCHEME"

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

            ./prep.sh \
                $SOLVER \
                $MESHX \
                $NPROCX \
                $PSOLVER \
                $NORMALSCHEME \
                $CURVATURESCHEME

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

            $PYTHON post.py log.$SOLVER

            if [ -f "result.txt" ]; then

                E1=$(sed -n '1p' < result.txt)
                E2=$(sed -n '2p' < result.txt)
                E3=$(sed -n '3p' < result.txt)
                E4=$(sed -n '4p' < result.txt)

                P1=$(sed -n '5p' < result.txt)
                P2=$(sed -n '6p' < result.txt)
                P3=$(sed -n '7p' < result.txt)
                P4=$(sed -n '8p' < result.txt)

                NDT=$(sed -n '9p' < result.txt)
                NITER=$(sed -n '10p' < result.txt)

            fi

            echo \
                "$SOLVER," \
                "$MESH," \
                "$NPROC," \
                "$PSOLVER," \
                "$NORMALSCHEME," \
                "$CURVATURESCHEME," \
                "$E1," \
                "$E2," \
                "$E3," \
                "$E4," \
                "$P1," \
                "$P2," \
                "$P3," \
                "$P4," \
                "$NDT," \
                "$NITER" \ >> $CURR/$CSV

            ##

            cd $CURR

            echo "Finished $CASE"

        ) &
    fi

done
done
done
done
done
done

wait

rm -f $TASKFILE.[0-9]+
