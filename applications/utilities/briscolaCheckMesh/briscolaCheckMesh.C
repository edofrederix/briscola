#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"

#include "fvMesh.H"
#include "meshFields.H"

#include "infoFunctions.H"
#include "checkFunctions.H"

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

    arguments::addBoolOption("brickInfo", "Show brick info");
    arguments::addBoolOption("patchInfo", "Show patch info");
    arguments::addBoolOption("levelInfo", "Show level info");
    arguments::addBoolOption("parallelInfo", "Show parallel info");
    arguments::addBoolOption("aggInfo", "Show level agglomerate info");

    arguments::addBoolOption("v", "Verbose info");

    arguments::addBoolOption
    (
        "parallelConnectivityInfo",
        "Show parallel connectivity info"
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

    const bool verbose = args.optionFound("v");

    generalInfo(fvMsh, verbose);

    if (args.optionFound("brickInfo"))
        brickInfo(fvMsh, verbose);

    if (args.optionFound("patchInfo"))
        patchInfo(fvMsh, verbose);

    if (args.optionFound("levelInfo"))
        levelInfo(fvMsh, verbose);

    if (Pstream::parRun() && args.optionFound("parallelInfo"))
        parallelInfo(fvMsh, verbose);

    if (Pstream::parRun() && args.optionFound("aggInfo"))
        aggInfo(fvMsh, verbose);

    if (Pstream::parRun())
    {
        if (args.optionFound("parallelConnectivityInfo"))
            parallelConnectivityInfo(fvMsh, verbose);

        returnReduce(0, sumOp<label>());

        checkParallelFaceCenters<colocated>(fvMsh);
        checkParallelFaceDeltas<colocated>(fvMsh);

        if (fvMsh.structured())
        {
            checkParallelFaceCenters<staggered>(fvMsh);
            checkParallelFaceDeltas<staggered>(fvMsh);
        }

        returnReduce(0, sumOp<label>());
    }

    Info<< nl << "All checks completed" << endl;
}
