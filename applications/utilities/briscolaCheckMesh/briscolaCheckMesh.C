#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fvMesh.H"
#include "meshFields.H"

#include "infoFunctions.H"
#include "parallelChecks.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    arguments::addOption
    (
        "dict",
        "file",
        "specify alternative dictionary for the briscolaMesh description"
    );

    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"

    fileName dictPath =
        args.optionFound("dict")
      ? args["dict"]
      : runTime.system()/"briscolaMeshDict";

    IOdictionary meshDict
    (
        IOobject
        (
            dictPath,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    fvMesh fvMsh(meshDict, runTime);

    generalInfo(fvMsh);
    brickInfo(fvMsh);
    patchInfo(fvMsh);
    parallelInfo(fvMsh);

    if (Pstream::parRun())
    {
        parallelConnectivityInfo(fvMsh);

        returnReduce(0, sumOp<label>());

        checkParallelFaceCenters<colocated>(fvMsh);
        checkParallelFaceDeltas<colocated>(fvMsh);

        // On triple or quintuple points parallel face and delta checks will
        // fail on staggered meshes

        if (fvMsh.topology().structured())
        {
            checkParallelFaceCenters<staggered>(fvMsh);
            checkParallelFaceDeltas<staggered>(fvMsh);
        }

        returnReduce(0, sumOp<label>());
    }
}
