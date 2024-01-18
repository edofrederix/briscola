#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"
#include "OSspecific.H"
#include "OFstream.H"

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

    arguments args(argc, argv);

    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    // Read arguments
    const label lineDir(args.argRead<label>(1));
    const word fieldName(args.argRead<word>(2));

    #include "createFields.H"
    #include "createBriscolaIO.H"

    if (lineDir < 0 || lineDir > 2)
    {
        FatalError << "Invalid line direction (0, 1, 2)" << endl;
        FatalError.exit();
    }

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
                avgVals[i] += field(0,0,i,j,k);
                break;

            case 1:
                avgVals[j] += field(0,0,i,j,k);
                break;

            case 2:
                avgVals[k] += field(0,0,i,j,k);
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
        forAll(coordinates, i)
        {
            coordinates[i] =
                fvMsh.msh().cast<rectilinearMesh>().globalPoints()[lineDir][i];
        }

        // Write line to file
        const fileName path("postProcessing/briscolaChannelPost"/fieldName);
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

        file << dir << "    " << fieldName << nl;

        forAll(globalAvgValsLine(), i)
        {
            file << coordinates[i] << "   ";
            file << globalAvgValsLine()[i] << nl;
        }
    }
}
