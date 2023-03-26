#include "rectilinearMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(rectilinearMesh, 0);
addToRunTimeSelectionTable(mesh, rectilinearMesh, dictionary);

void rectilinearMesh::setCellSizes()
{
    cellSizes_.clear();
    cellSizes_.setSize(3);

    for (int dir = 0; dir < 3; dir++)
    {
        const decompositionMap& map = decomp().map();

        cellSizes_.set(dir, new scalarList(N()[dir], 0.0));

        scalarList& sizes = cellSizes_[dir];

        labelVector base = units[dir];
        label cursor = 0;

        for (int i = 0; i < map.shape()[dir]; i++)
        {
            const label proc = map(base*i);

            if (proc == Pstream::myProcNo())
            {
                scalarList partSizes
                (
                    this->operator[](0).rectilinearCellSizes(dir)
                );

                if (proc == Pstream::masterNo())
                {
                    // Copy directly

                    forAll(partSizes, j)
                    {
                        sizes[cursor+j] = partSizes[j];
                    }

                    cursor += partSizes.size();
                }
                else
                {
                    // Send to master

                    OPstream send
                    (
                        Pstream::commsTypes::blocking,
                        Pstream::masterNo()
                    );

                    send << partSizes;
                }
            }
            else if (Pstream::master())
            {
                // Receive from slave

                scalarList partSizes;

                IPstream recv
                (
                    Pstream::commsTypes::blocking,
                    proc
                );

                recv >> partSizes;

                forAll(partSizes, j)
                {
                    sizes[cursor+j] = partSizes[j];
                }

                cursor += partSizes.size();
            }
        }

        // Send back to slaves

        if (Pstream::master())
        {
            for (int proc = 0; proc < Pstream::nProcs(); proc++)
            if (proc != Pstream::masterNo())
            {
                OPstream send
                (
                    Pstream::commsTypes::blocking,
                    proc
                );

                send << sizes;
            }
        }
        else
        {
            IPstream recv
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            recv >> sizes;
        }
    }
}

rectilinearMesh::rectilinearMesh(const IOdictionary& dict)
:
    structuredMesh(dict)
{
    setCellSizes();
}

rectilinearMesh::rectilinearMesh(autoPtr<mesh>& mshPtr)
:
    structuredMesh(mshPtr)
{
    setCellSizes();
}

rectilinearMesh::rectilinearMesh(const rectilinearMesh& msh)
:
    structuredMesh(msh),
    cellSizes_(msh.cellSizes_)
{}

rectilinearMesh::rectilinearMesh(rectilinearMesh& msh, bool reuse)
:
    structuredMesh(msh, reuse),
    cellSizes_(msh.cellSizes_, reuse)
{}

rectilinearMesh::~rectilinearMesh()
{}

}

}
