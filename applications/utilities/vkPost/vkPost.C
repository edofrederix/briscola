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
#include "meshField.H"
#include "meshType.H"
#include "pointDataExchange.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");
    arguments::validArgs.append("x-coordinate");
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
    const scalar xCoord(args.argRead<scalar>(1));
    const word fieldName(args.argRead<word>(2));
    const bool stag(args.optionFound("staggered"));

    forAll(times, timei)
    {
        runTime.setTime(times[timei], timei);

        Info<< "Running vkPost for time = "
            << runTime.timeName() << endl;

        #include "createFields.H"
        #include "createBriscolaIO.H"

        if (xCoord < -15 || xCoord > 25)
        {
            FatalError << "Invalid x-coordinate, range=(-15, 25)" << endl;
            FatalError.exit();
        }

        if (!stag)
        {
            // Mesh dimensions of processor part
            labelVector N(fvMsh.N<colocated>());

            // Global mesh dimensions
            labelVector Ng(fvMsh.msh().cast<rectilinearMesh>().N());

            // Averaging direction coordinates
            scalarList zCoordinates(Ng.z(), Zero);

            const PartialList<scalar>& zGlobalPoints =
                fvMsh.msh().cast<rectilinearMesh>().globalPoints()[2];

            forAll(zCoordinates, i)
            {
                zCoordinates[i] =
                    (
                        zGlobalPoints[i]
                        + zGlobalPoints[i+1]
                    )
                    /
                    2.0;
            }

            // Averaged line coordinates
            scalarList yCoordinates(Ng.y(), Zero);

            const PartialList<scalar>& yGlobalPoints =
                fvMsh.msh().cast<rectilinearMesh>().globalPoints()[1];

            forAll(yCoordinates, i)
            {
                yCoordinates[i] =
                    (
                        yGlobalPoints[i]
                        + yGlobalPoints[i+1]
                    )
                    /
                    2.0;
            }

            // Moving averaging line coordinates
            vectorList coordinates
            (
                Pstream::master() ? Ng.y() : 0,
                Zero
            );

            forAll(coordinates, i)
            {
                coordinates[i] = vector
                (
                    xCoord,
                    yCoordinates[i],
                    zCoordinates[0]
                );
            }

            // List of averaged values
            scalarList colAveragedLine(Ng.y(), Zero);

            forAll(zCoordinates, i)
            {
                forAll(coordinates, j)
                {
                    coordinates[j].z() = zCoordinates[i];
                }

                pointDataExchange<colocated> exchange
                (
                    coordinates, fvMsh, 0, 0
                );

                scalarList exchangeData(move(exchange(colField)));

                forAll(exchangeData, j)
                {
                    colAveragedLine[j] += exchangeData[j];
                }
            }

            forAll(colAveragedLine, i)
            {
                colAveragedLine[i] /= Ng.z();
            }

            if (Pstream::master())
            {
                // Write line to file
                const fileName path
                (
                    "postProcessing/vkPost"
                    /runTime.timeName()
                    /colField.name()
                    /name(xCoord)
                );
                mkDir(path);

                OFstream file(path/"averagedLine.txt");

                file << "y;" << colField.name() << nl;

                forAll(colAveragedLine, i)
                {
                    file << yCoordinates[i] << ";";
                    file << colAveragedLine[i] << nl;
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
                if (d == 1)
                {
                    Ng[d] += 1;
                }

                // Averaging direction coordinates
                scalarList zCoordinates(Ng.z(), Zero);

                const PartialList<scalar>& zGlobalPoints =
                    fvMsh.msh().cast<rectilinearMesh>().globalPoints()[2];

                forAll(zCoordinates, i)
                {
                    zCoordinates[i] =
                        (
                            zGlobalPoints[i]
                            + zGlobalPoints[i+1]
                        )
                        /
                        2.0;
                }

                // Averaged line coordinates
                scalarList yCoordinates(Ng.y(), Zero);

                const PartialList<scalar>& yGlobalPoints =
                    fvMsh.msh().cast<rectilinearMesh>().globalPoints()[1];

                forAll(yCoordinates, i)
                {
                    if (d == 1)
                    {
                        yCoordinates[i] = yGlobalPoints[i];
                    }
                    else
                    {
                        yCoordinates[i] =
                            (
                                yGlobalPoints[i]
                                + yGlobalPoints[i+1]
                            )
                            /
                            2.0;
                    }
                }

                // Moving averaging line coordinates
                vectorList coordinates
                (
                    Pstream::master() ? Ng.y() : 0,
                    Zero
                );

                forAll(coordinates, i)
                {
                    coordinates[i] = vector
                    (
                        xCoord,
                        yCoordinates[i],
                        zCoordinates[0]
                    );
                }

                // List of averaged values
                scalarList stagAveragedLine(Ng.y(), Zero);

                forAll(zCoordinates, i)
                {
                    forAll(coordinates, j)
                    {
                        coordinates[j].z() = zCoordinates[i];
                    }

                    pointDataExchange<staggered> exchange
                    (
                        coordinates, fvMsh, 0, d
                    );

                    scalarList exchangeData(move(exchange(stagField)));

                    forAll(exchangeData, j)
                    {
                        stagAveragedLine[j] += exchangeData[j];
                    }
                }

                forAll(stagAveragedLine, i)
                {
                    stagAveragedLine[i] /= Ng.z();
                }

                if (Pstream::master())
                {
                    // Write line to file
                    const fileName path
                    (
                        "postProcessing/vkPost"
                        /runTime.timeName()
                        /stagField.name()+"_"+Foam::name(d)
                        /name(xCoord)
                    );
                    mkDir(path);

                    OFstream file(path/"averagedLine.txt");

                    file << "y;" << stagField.name()+"_"+Foam::name(d) << nl;

                    forAll(stagAveragedLine, i)
                    {
                        file << yCoordinates[i] << ";";
                        file << stagAveragedLine[i] << nl;
                    }
                }
            }
        }
    }
}
