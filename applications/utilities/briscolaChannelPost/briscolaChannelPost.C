#include "arguments.H"
#include "IOdictionary.H"
#include "Time.H"
#include "timeSelector.H"
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

template<class Type>
void printComponents(OFstream& file, const Type& val)
{
    NotImplemented::error("printComponents", "Unsupported type");
}

template<>
void printComponents<scalar>(OFstream& file, const scalar& val)
{
    file << val << nl;
}

template<>
void printComponents<vector>(OFstream& file, const vector& val)
{
    file << val.x() << ";" << val.y() << ";" << val.z() << nl;
}

template<>
void printComponents<symmTensor>(OFstream& file, const symmTensor& val)
{
    file << val.xx() << ";" << val.xy() << ";" << val.xz() << ";"
         << val.yy() << ";" << val.yz() << ";" << val.zz() << nl;
}

template<>
void printComponents<tensor>(OFstream& file, const tensor& val)
{
    file << val.xx() << ";" << val.xy() << ";" << val.xz() << ";"
         << val.yx() << ";" << val.yy() << ";" << val.yz() << ";"
         << val.zx() << ";" << val.zy() << ";" << val.zz() << nl;
}

template<class Type, class MeshType>
int briscolaChannelPost(arguments args)
{
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

    forAll(times, timei)
    {
        runTime.setTime(times[timei], timei);

        Info<< "Running briscolaChannelPost for time = "
            << runTime.name() << endl;

        meshField<Type,MeshType> field
        (
            fieldName,
            fvMsh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            true
        );

        #include "createBriscolaIO.H"

        if (lineDir < 0 || lineDir > 2)
        {
            FatalError << "Invalid line direction (0, 1, 2)" << endl;
            FatalError.exit();
        }

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            // Mesh dimensions of processor part
            labelVector N(fvMsh.N<MeshType>(d));

            // Last processor index
            labelVector L(fvMsh.msh().decomp().map().legend()[Pstream::nProcs() - 1]);

            // Global mesh dimensions
            labelVector Ng(fvMsh.msh().cast<rectilinearMesh>().N());
            if (word(MeshType::typeName) == "staggered" && d == lineDir)
            {
                Ng[d] += 1;
            }

            // Local list of averaged (summed) values
            List<Type> avgVals(N[lineDir], Zero);

            // Sum up the values on each processor
            for (int i = 0; i < N.x(); i++)
            for (int j = 0; j < N.y(); j++)
            for (int k = 0; k < N.z(); k++)
            {
                switch (lineDir)
                {
                    case 0:
                        avgVals[i] += field(0,d,i,j,k);
                        break;

                    case 1:
                        avgVals[j] += field(0,d,i,j,k);
                        break;

                    case 2:
                        avgVals[k] += field(0,d,i,j,k);
                        break;
                }
            }

            // Contains all local summed lines one after the other
            autoPtr<List<Type>> globalAvgValsList;
            // Global spatially averaged line
            autoPtr<List<Type>> globalAvgValsLine;

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
                           (word(MeshType::typeName) == "staggered")
                        && (d == lineDir)
                        && (fvMsh.msh().decomp().map().legend()[p][lineDir] == L[lineDir])
                    )
                    {
                        globalListSize++;
                    }
                }

                globalAvgValsList.reset
                (
                    new List<Type>
                    (
                        globalListSize,
                        Zero
                    )
                );
                globalAvgValsLine.reset
                (
                    new List<Type>
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
                    .partSizePerProc()[p][lineDir]*sizeof(Type);

                if
                (
                       (word(MeshType::typeName) == "staggered")
                    && (d == lineDir)
                    && (fvMsh.msh().decomp().map().legend()[p][lineDir] == L[lineDir])
                )
                {
                    rs[p] += sizeof(Type);
                }

                displacement += rs[p];
            }

            // Gather averaged (summed) lists on master processor
            UPstream::gather
            (
                reinterpret_cast<char*>(avgVals.begin()),
                avgVals.size()*sizeof(Type),
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
                    rd[p] /= sizeof(Type);
                    rs[p] /= sizeof(Type);
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
                    if (word(MeshType::typeName) == "staggered" && d == lineDir)
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

                word suffix
                (
                    word(MeshType::typeName) == "staggered" ?
                    "_"+Foam::name(d)
                    : ""
                );

                // Write line to file
                const fileName path
                (
                    "postProcessing/briscolaChannelPost"
                    /runTime.name()
                    /fieldName+suffix
                );
                mkDir(path);

                OFstream file(path/"averagedLine.txt");

                forAll(globalAvgValsLine(), i)
                {
                    file << coordinates[i] << ";";
                    printComponents(file, globalAvgValsLine()[i]);
                }
            }
        }
    }

    return 0;
}

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");
    arguments::validArgs.append("line direction");
    arguments::validArgs.append("field to average");
    arguments::validArgs.append("field type");
    arguments::addOption
    (
        "time",
        "ranges",
        "comma-separated time ranges - eg, ':10,20,40:70,1000:'"
    );

    arguments args(argc, argv);

    const word type(args.argRead<word>(3));

    if (type == "colocatedScalar")
        return briscolaChannelPost<scalar,colocated>(args);
    else if (type == "colocatedVector")
        return briscolaChannelPost<vector,colocated>(args);
    else if (type == "colocatedTensor")
        return briscolaChannelPost<tensor,colocated>(args);
    else if (type == "colocatedSymmTensor")
        return briscolaChannelPost<symmTensor,colocated>(args);
    else if (type == "staggeredScalar")
        return briscolaChannelPost<scalar,staggered>(args);
    else
    {
        FatalError
                << "Type " << type
                << " not supported."
                << endl << abort(FatalError);

        return 1;
    }
}
