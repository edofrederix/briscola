#include "rectilinearMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(rectilinearMesh, 0);
addToRunTimeSelectionTable(mesh, rectilinearMesh, dictionary);

void rectilinearMesh::setMetrics()
{
    const decompositionMap& map = decomp().map();
    const partLevel& part = this->operator[](0);
    const partLevelPoints& points = part.points();

    // Base vectors

    vectorList base(3);

    for (int d = 0; d < 3; d++)
    {
        base[d] = points(units[d]) - points(zeroXYZ);
        base[d] /= Foam::mag(base[d]);
    }

    base_ = tensor(base[0], base[1], base[2]);

    if (base_.x() != base[0] || base_.y() != base[1] || base_.z() != base[2])
    {
        FatalErrorInFunction
            << "Base tensor does not match base vectors" << endl
            << abort(FatalError);
    }

    // Local cell sizes and points

    localCellSizes_.clear();
    localPoints_.clear();
    localCellSizes_.setSize(3);
    localPoints_.setSize(3);

    for (int d = 0; d < 3; d++)
    {
        localCellSizes_.set(d, new scalarList(part.N()[d]));
        localPoints_.set(d, new scalarList(part.N()[d]+1));

        scalarList& localCellSizes = localCellSizes_[d];
        scalarList& localPoints = localPoints_[d];

        const labelVector dir = units[d];

        forAll(localCellSizes, i)
        {
            localCellSizes[i] =
                Foam::mag(points(dir*(i+1)) - points(dir*i));
        }

        const vector base =
            d == 0 ? base_.x() : d == 1 ? base_.y() : base_.z();

        forAll(localPoints, i)
        {
            localPoints[i] = points(dir*i) & base;
        }
    }

    // Global cell sizes

    globalCellSizes_.clear();
    globalPoints_.clear();
    globalCellSizes_.setSize(3);
    globalPoints_.setSize(3);

    for (int d = 0; d < 3; d++)
    {
        globalCellSizes_.set(d, new scalarList(this->N()[d]));
        globalPoints_.set(d, new scalarList(this->N()[d]+1));

        scalarList& globalCellSizes = globalCellSizes_[d];
        scalarList& globalPoints = globalPoints_[d];

        const labelVector base = units[d];
        label cursor = 0;

        for (int i = 0; i < map.shape()[d]; i++)
        {
            const label proc = map(base*i);

            if (proc == Pstream::myProcNo())
            {
                const scalarList& localCellSizes = localCellSizes_[d];
                const scalarList& localPoints = localPoints_[d];

                if (proc == Pstream::masterNo())
                {
                    // Copy directly

                    forAll(localCellSizes, j)
                    {
                        globalCellSizes[cursor+j] = localCellSizes[j];
                    }

                    forAll(localPoints, j)
                    {
                        globalPoints[cursor+j] = localPoints[j];
                    }

                    cursor += localCellSizes.size();
                }
                else
                {
                    // Send to master

                    OPstream send
                    (
                        Pstream::commsTypes::blocking,
                        Pstream::masterNo()
                    );

                    send << localCellSizes;
                    send << localPoints;
                }
            }
            else if (Pstream::master())
            {
                // Receive from slave

                scalarList localCellSizes;
                scalarList localPoints;

                IPstream recv
                (
                    Pstream::commsTypes::blocking,
                    proc
                );

                recv >> localCellSizes;
                recv >> localPoints;

                forAll(localCellSizes, j)
                {
                    globalCellSizes[cursor+j] = localCellSizes[j];
                }

                forAll(localPoints, j)
                {
                    globalPoints[cursor+j] = localPoints[j];
                }

                cursor += localCellSizes.size();
            }
        }

        // Send back to slaves

        Pstream::scatter(globalCellSizes);
        Pstream::scatter(globalPoints);
    }
}

rectilinearMesh::rectilinearMesh(const IOdictionary& dict)
:
    structuredMesh(dict)
{
    setMetrics();
}

rectilinearMesh::rectilinearMesh(autoPtr<mesh>& mshPtr)
:
    structuredMesh(mshPtr)
{
    setMetrics();
}

rectilinearMesh::rectilinearMesh(const rectilinearMesh& msh)
:
    structuredMesh(msh),
    base_(msh.base_),
    localCellSizes_(msh.localCellSizes_),
    localPoints_(msh.localPoints_),
    globalCellSizes_(msh.globalCellSizes_),
    globalPoints_(msh.globalPoints_)
{}

rectilinearMesh::rectilinearMesh(rectilinearMesh& msh, bool reuse)
:
    structuredMesh(msh, reuse),
    base_(msh.base_),
    localCellSizes_(msh.localCellSizes_),
    localPoints_(msh.localPoints_),
    globalCellSizes_(msh.globalCellSizes_),
    globalPoints_(msh.globalPoints_)
{}

rectilinearMesh::~rectilinearMesh()
{}

}

}
