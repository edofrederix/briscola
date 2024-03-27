#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"
#include "timeSelector.H"
#include "OSspecific.H"
#include "OFstream.H"
#include "PtrList.H"

#include "fv.H"
#include "fvMesh.H"
#include "rectilinearMesh.H"
#include "meshFields.H"
#include "meshType.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");
    arguments::validArgs.append("line direction");
    arguments::validArgs.append("field to average");
    arguments::addBoolOption("staggered", "run for staggered field");
    arguments::addOption
    (
        "time",
        "ranges",
        "comma-separated time ranges - eg, ':10,20,40:70,1000:'"
    );

    arguments args(argc, argv);

    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    Foam::instantList timeDirs = runTime.times();
    List<bool> selectTimes(timeDirs.size(), true);

    selectTimes = Foam::timeSelector
    (
        args.optionLookup("time")()
    ).selected(timeDirs);

    Foam::instantList times(subset(selectTimes, timeDirs));

    // Read arguments
    const label lineDir(args.argRead<label>(1));
    const word fieldName(args.argRead<word>(2));
    const bool stag(args.optionFound("staggered"));

    forAll(times, timei)
    {
        runTime.setTime(times[timei], timei);

        Info<< "Running briscolaChannelPost for time = "
            << runTime.timeName() << endl;

        #include "createFields.H"
        #include "createBriscolaIO.H"

        if (lineDir < 0 || lineDir > 2)
        {
            FatalError << "Invalid line direction (0, 1, 2)" << endl;
            FatalError.exit();
        }

        if (!stag)
        {
            // Mesh dimensions of processor part
            labelVector N(fvMsh.N<colocated>());

            // Global mesh dimensions
            labelVector Ng(fvMsh.msh().cast<rectilinearMesh>().N());

            // Local list of averaged (summed) values
            scalarList avgVals(N[lineDir], Zero);

            // Sum up the values on each processor
            for (int i = 0; i < N.x(); i++)
            for (int j = 0; j < N.y(); j++)
            for (int k = 0; k < N.z(); k++)
            {
                switch (lineDir)
                {
                    case 0:
                        avgVals[i] += colField(0,0,i,j,k);
                        break;

                    case 1:
                        avgVals[j] += colField(0,0,i,j,k);
                        break;

                    case 2:
                        avgVals[k] += colField(0,0,i,j,k);
                        break;
                }
            }

            // Contains all local summed lines one after the other
            autoPtr<scalarList> globalAvgValsList;
            // Global spatially averaged line
            autoPtr<scalarList> globalAvgValsLine;

            // Global lists only needed on master proc
            if (Pstream::master())
            {
                // Set global list size
                label globalListSize = 0;
                for (int p = 0; p < Pstream::nProcs(); p++)
                {
                    globalListSize +=
                        fvMsh.msh().decomp().partSizePerProc()[p][lineDir];
                }

                globalAvgValsList.reset(new scalarList(globalListSize, Zero));
                globalAvgValsLine.reset
                (
                    new scalarList
                    (
                        Ng[lineDir],
                        Zero
                    )
                );
            }

            // Receive sizes and displacements
            labelList rs(Pstream::nProcs(), Zero);
            labelList rd(Pstream::nProcs(), Zero);

            label displacement = 0;
            forAll(rs, p)
            {
                rd[p] = displacement;
                rs[p] = fvMsh.msh().decomp()
                    .partSizePerProc()[p][lineDir]*sizeof(scalar);
                displacement += rs[p];
            }

            // Gather averaged (summed) lists on master processor
            UPstream::gather
            (
                reinterpret_cast<char*>(avgVals.begin()),
                avgVals.size()*sizeof(scalar),
                reinterpret_cast<char*>
                (
                    Pstream::master()
                ? globalAvgValsList->begin()
                : nullptr
                ),
                rs,
                rd,
                UPstream::worldComm
            );

            if (Pstream::master())
            {
                // Reset rd and rs from bytes to number of elements
                forAll(rs, p)
                {
                    rd[p] /= sizeof(scalar);
                    rs[p] /= sizeof(scalar);
                }

                for (int p = 0; p < Pstream::nProcs(); p++)
                {
                    // processor p starting index
                    label startIndex =
                        fvMsh.msh().decomp().globalStartPerProc()[p][lineDir];

                    for (int i = 0; i < rs[p]; i++)
                    {
                        globalAvgValsLine()[startIndex+i]
                            += globalAvgValsList()[rd[p]+i];
                    }
                }

                label avgN = cmptProduct(Ng);
                avgN /= Ng[lineDir];

                forAll(globalAvgValsLine(), i)
                {
                    globalAvgValsLine()[i] /= avgN;
                }

                // Line point coordinates
                scalarList coordinates(Ng[lineDir], Zero);
                const PartialList<scalar>& globalPoints =
                    fvMsh.msh().cast<rectilinearMesh>().globalPoints()[lineDir];

                forAll(coordinates, i)
                {
                    coordinates[i] =
                        (
                            globalPoints[i]
                            + globalPoints[i+1]
                        )
                        /
                        2.0;
                }

                // Write line to file
                const fileName path
                (
                    "postProcessing/briscolaChannelPost"
                    /runTime.timeName()
                    /colField.name()
                );
                mkDir(path);

                OFstream file(path/"averagedLine.txt");

                word dir(" ");
                switch (lineDir)
                {
                    case 0:
                        dir = "x";
                        break;

                    case 1:
                        dir = "y";
                        break;

                    case 2:
                        dir = "z";
                        break;
                }

                file << dir << ";" << colField.name() << nl;

                forAll(globalAvgValsLine(), i)
                {
                    file << coordinates[i] << ";";
                    file << globalAvgValsLine()[i] << nl;
                }
            }
        }
        else
        {
            for (int d = 0; d < 3; d++)
            {
                // Mesh dimensions of processor part
                labelVector N(fvMsh.N<staggered>(d));

                // Last processor index
                labelVector L(fvMsh.msh().decomp().map().legend()[Pstream::nProcs() - 1]);

                // Global mesh dimensions
                labelVector Ng(fvMsh.msh().cast<rectilinearMesh>().N());
                if (d == lineDir)
                {
                    Ng[d] += 1;
                }

                // Local list of averaged (summed) values
                scalarList avgVals(N[lineDir], Zero);

                // Sum up the values on each processor
                for (int i = 0; i < N.x(); i++)
                for (int j = 0; j < N.y(); j++)
                for (int k = 0; k < N.z(); k++)
                {
                    switch (lineDir)
                    {
                        case 0:
                            avgVals[i] += stagField(0,d,i,j,k);
                            break;

                        case 1:
                            avgVals[j] += stagField(0,d,i,j,k);
                            break;

                        case 2:
                            avgVals[k] += stagField(0,d,i,j,k);
                            break;
                    }
                }

                // Contains all local summed lines one after the other
                autoPtr<scalarList> globalAvgValsList;
                // Global spatially averaged line
                autoPtr<scalarList> globalAvgValsLine;

                // Global lists only needed on master proc
                if (Pstream::master())
                {
                    // Set global list size
                    label globalListSize = 0;
                    for (int p = 0; p < Pstream::nProcs(); p++)
                    {
                        globalListSize +=
                            fvMsh.msh().decomp().partSizePerProc()[p][lineDir];

                        if
                        (
                            (d == lineDir)
                            && (fvMsh.msh().decomp().map().legend()[p][lineDir] == L[lineDir])
                        )
                        {
                            globalListSize++;
                        }
                    }

                    globalAvgValsList.reset
                    (
                        new scalarList
                        (
                            globalListSize,
                            Zero
                        )
                    );
                    globalAvgValsLine.reset
                    (
                        new scalarList
                        (
                            Ng[lineDir],
                            Zero
                        )
                    );
                }

                // Receive sizes and displacements
                labelList rs(Pstream::nProcs(), Zero);
                labelList rd(Pstream::nProcs(), Zero);

                label displacement = 0;
                forAll(rs, p)
                {
                    rd[p] = displacement;
                    rs[p] = fvMsh.msh().decomp()
                        .partSizePerProc()[p][lineDir]*sizeof(scalar);

                    if
                    (
                        (d == lineDir)
                        && (fvMsh.msh().decomp().map().legend()[p][lineDir] == L[lineDir])
                    )
                    {
                        rs[p] += sizeof(scalar);
                    }

                    displacement += rs[p];
                }

                // Gather averaged (summed) lists on master processor
                UPstream::gather
                (
                    reinterpret_cast<char*>(avgVals.begin()),
                    avgVals.size()*sizeof(scalar),
                    reinterpret_cast<char*>
                    (
                        Pstream::master()
                    ? globalAvgValsList->begin()
                    : nullptr
                    ),
                    rs,
                    rd,
                    UPstream::worldComm
                );

                if (Pstream::master())
                {
                    // Reset rd and rs from bytes to number of elements
                    forAll(rs, p)
                    {
                        rd[p] /= sizeof(scalar);
                        rs[p] /= sizeof(scalar);
                    }

                    for (int p = 0; p < Pstream::nProcs(); p++)
                    {
                        // processor p starting index
                        label startIndex =
                            fvMsh.msh().decomp().globalStartPerProc()[p][lineDir];

                        scalar s = rs[p];

                        for (int i = 0; i < s; i++)
                        {
                            globalAvgValsLine()[startIndex+i]
                                += globalAvgValsList()[rd[p]+i];
                        }
                    }

                    label avgN = cmptProduct(Ng);
                    avgN /= Ng[lineDir];

                    forAll(globalAvgValsLine(), i)
                    {
                        globalAvgValsLine()[i] /= avgN;
                    }

                    // Line point coordinates
                    scalarList coordinates(Ng[lineDir], Zero);
                    const PartialList<scalar>& globalPoints =
                        fvMsh.msh().cast<rectilinearMesh>().globalPoints()[lineDir];

                    forAll(coordinates, i)
                    {
                        if (d == lineDir)
                        {
                            coordinates[i] = globalPoints[i];
                        }
                        else
                        {
                            coordinates[i] =
                            (
                                globalPoints[i]
                                + globalPoints[i+1]
                            )
                            /
                            2.0;
                        }
                    }

                    // Write line to file
                    const fileName path
                    (
                        "postProcessing/briscolaChannelPost"
                        /runTime.timeName()
                        /stagField.name()+"_"+Foam::name(d)
                    );
                    mkDir(path);

                    OFstream file(path/"averagedLine.txt");

                    word dir(" ");
                    switch (lineDir)
                    {
                        case 0:
                            dir = "x";
                            break;

                        case 1:
                            dir = "y";
                            break;

                        case 2:
                            dir = "z";
                            break;
                    }

                    file << dir << ";" << stagField.name()+"_"+Foam::name(d) << nl;

                    forAll(globalAvgValsLine(), i)
                    {
                        file << coordinates[i] << ";";
                        file << globalAvgValsLine()[i] << nl;
                    }
                }
            }
        }
    }
}
