#!/bin/bash

# By default 64 cores are used. Simulations are sequentially started only when
# resources are available. Output is written to results.csv.

NTASKS="${SLURM_NTASKS:-64}"

# Some parameters

RUNDIR="runs"
CSV="results.csv"
TEMPLATE="template"
PYTHON="python3"

# Simulation parameters

SOLVERS=(briscolaColocatedTwoPhase briscolaStaggeredTwoPhase)
MESHES=(64 128)
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
##

source ../procManagement.sh

for I in "${!SOLVERS[@]}"; do
for J in "${!MESHES[@]}"; do
for K in "${!PSOLVERS[@]}"; do
for L in "${!NORMALSCHEMES[@]}"; do
for M in "${!CURVATURESCHEMES[@]}"; do

    sleep 1

    SOLVER=${SOLVERS[$I]}
    MESH=${MESHES[$J]}
    NPROCPERBRICKSIDE=${NPROCSPERBRICKSIDE[$J]}
    PSOLVER=${PSOLVERS[$K]}
    NORMALSCHEME=${NORMALSCHEMES[$L]}
    CURVATURESCHEME=${CURVATURESCHEMES[$M]}

    NPROCX=$NPROCPERBRICKSIDE

    NPROC=$(echo "print(int($NPROCX**3))" | python)

    MESHX=$MESH

    CASE="$SOLVER-$MESH-$NPROC-$PSOLVER-$NORMALSCHEME-$CURVATURESCHEME"

    wait_for_procs $NPROC $NTASKS

    (
        echo "Starting $CASE"

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
            srun --exclusive -n $NPROC $SOLVER > log.$SOLVER
        else
            srun --exclusive -n $NPROC $SOLVER -parallel > log.$SOLVER
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

    store_procs $NPROC $!

done
done
done
done
done

wait
