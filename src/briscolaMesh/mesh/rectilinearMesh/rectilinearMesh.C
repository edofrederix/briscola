#include "rectilinearMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "parallelBoundary.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(rectilinearMesh, 0);
addToRunTimeSelectionTable(mesh, rectilinearMesh, dictionary);

void rectilinearMesh::setMetrics()
{
    const decompositionMap& map = decomp().map();
    const part& p = this->operator[](0);
    const partPoints& points = p.points();

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

    localCellSizesData_.clear();
    localPointsData_.clear();
    localCellSizesData_.setSize(3);
    localPointsData_.setSize(3);

    for (int d = 0; d < 3; d++)
    {
        localCellSizesData_.set(d, new scalarList(p.N()[d]+2));
        localPointsData_.set(d, new scalarList(p.N()[d]+3));

        scalarList& localCellSizesData = localCellSizesData_[d];
        scalarList& localPointsData = localPointsData_[d];

        const labelVector dir = units[d];

        const vector base =
            d == 0 ? base_.x() : d == 1 ? base_.y() : base_.z();

        forAll(localPointsData, i)
        {
            localPointsData[i] = trimPrecision(points(dir*(i-1)) & base);
        }

        forAll(localCellSizesData, i)
        {
            localCellSizesData[i] =
                localPointsData[i+1] - localPointsData[i];
        }

        // Check if local point coordinates form a monotonically increasing list

        for (label i = 1; i < localPointsData.size(); i++)
        {
            if (localPointsData[i] <= localPointsData[i-1])
            {
                FatalErrorInFunction
                    << "Could not generate local point list on rectilinear mesh"
                    << endl << abort(FatalError);
            }
        }
    }

    // Create local lists without ghosts as sub-lists

    localCellSizes_.clear();
    localPoints_.clear();
    localCellSizes_.setSize(3);
    localPoints_.setSize(3);

    for (int d = 0; d < 3; d++)
    {
        localCellSizes_.set
        (
            d,
            new PartialList<scalar>
            (
                localCellSizesData_[d],
                localCellSizesData_[d].size()-2,
                1
            )
        );

        localPoints_.set
        (
            d,
            new PartialList<scalar>
            (
                localPointsData_[d],
                localPointsData_[d].size()-2,
                1
            )
        );
    }

    // Global cell sizes

    globalCellSizesData_.clear();
    globalPointsData_.clear();
    globalCellSizesData_.setSize(3);
    globalPointsData_.setSize(3);

    PtrList<labelList> starts(3);

    for (int d = 0; d < 3; d++)
    {
        globalCellSizesData_.set(d, new scalarList(this->N()[d]+2));
        globalPointsData_.set(d, new scalarList(this->N()[d]+3));

        scalarList& globalCellSizes = globalCellSizesData_[d];
        scalarList& globalPoints = globalPointsData_[d];

        const labelVector base = units[d];

        label cursor = 0;

        starts.set(d, new labelList(map.shape()[d]));

        for (int i = 0; i < map.shape()[d]; i++)
        {
            const label proc = map(base*i);

            if (proc == Pstream::myProcNo())
            {
                const scalarList& localCellSizesData =
                    localCellSizesData_[d];

                const scalarList& localPointsData =
                    localPointsData_[d];

                if (proc == Pstream::masterNo())
                {
                    // Copy directly

                    forAll(localCellSizesData, j)
                    {
                        globalCellSizes[cursor+j] = localCellSizesData[j];
                    }

                    forAll(localPointsData, j)
                    {
                        globalPoints[cursor+j] = localPointsData[j];
                    }

                    starts[d][i] = cursor;
                    cursor += localCellSizesData.size()-2;
                }
                else
                {
                    // Send to master

                    OPstream send
                    (
                        Pstream::commsTypes::blocking,
                        Pstream::masterNo()
                    );

                    send << localCellSizesData;
                    send << localPointsData;
                }
            }
            else if (Pstream::master())
            {
                // Receive from slave

                scalarList localCellSizesData;
                scalarList localPointsData;

                IPstream recv
                (
                    Pstream::commsTypes::blocking,
                    proc
                );

                recv >> localCellSizesData;
                recv >> localPointsData;

                forAll(localCellSizesData, j)
                {
                    globalCellSizes[cursor+j] = localCellSizesData[j];
                }

                forAll(localPointsData, j)
                {
                    globalPoints[cursor+j] = localPointsData[j];
                }

                starts[d][i] = cursor;
                cursor += localCellSizesData.size()-2;
            }
        }

        // Send back to slaves

        Pstream::scatter(globalCellSizes);
        Pstream::scatter(globalPoints);

        if (Pstream::master())
        {
            for (int proc = 0; proc < Pstream::nProcs(); proc++)
            if (proc != Pstream::myProcNo())
            {
                OPstream send
                (
                    Pstream::commsTypes::blocking,
                    proc
                );

                send << starts[d];
            }
        }
        else
        {
            IPstream recv
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            recv >> starts[d];
        }
    }

    // Create global lists without ghosts as sub-lists

    globalCellSizes_.clear();
    globalPoints_.clear();
    globalCellSizes_.setSize(3);
    globalPoints_.setSize(3);

    for (int d = 0; d < 3; d++)
    {
        globalCellSizes_.set
        (
            d,
            new PartialList<scalar>
            (
                globalCellSizesData_[d],
                globalCellSizesData_[d].size()-2,
                1
            )
        );

        globalPoints_.set
        (
            d,
            new PartialList<scalar>
            (
                globalPointsData_[d],
                globalPointsData_[d].size()-2,
                1
            )
        );
    }

    // Set global starts

    globalStarts_.setSize(Pstream::nProcs());

    forAllBlock(map, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const label proc = map(ijk);

        globalStarts_[proc] = Zero;

        for (int d = 0; d < 3; d++)
            globalStarts_[proc][d] = starts[d][ijk[d]];
    }

    // Check if global mesh directions are uniform

    globalUniform_ = labelVector(1,1,1);

    const scalar tol = 1e-7;

    forAll(globalCellSizes_, d)
    {
        scalar maxCellSize = max(globalCellSizes_[d]);
        scalar minCellSize = min(globalCellSizes_[d]);

        if (maxCellSize - minCellSize > tol)
        {
            globalUniform_[d] = 0;
        }
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
    globalPoints_(msh.globalPoints_),
    globalUniform_(msh.globalUniform_)
{}

rectilinearMesh::~rectilinearMesh()
{}

labelVector rectilinearMesh::findCell(const vector& p, const label l) const
{
    const vector q
    (
        trimPrecision
        (
            vector
            (
                p & base_.x(),
                p & base_.y(),
                p & base_.z()
            )
        )
    );

    const PtrList<PartialList<scalar>>& x = localPoints_;

    // Check if the point is outside the domain boundaries

    for (label d = 0; d < 3; d++)
        if (q[d] < x[d].first() || q[d] > x[d].last() + 1e-12)
            return -unitXYZ;

    // Check if the point is just above/below the upper/lower boundary. If so,
    // the point is included if the boundary is not be a parallel one. If the
    // point is exactly on the lower boundary it must be included anyway.

    for (label d = 0; d < 3; d++)
        if (q[d] >= x[d].last() && q[d] <= x[d].last() + 1e-12)
            if (b(units[d]).castable<parallelBoundary>())
                return -unitXYZ;

    for (label d = 0; d < 3; d++)
        if (q[d] < x[d].first() && q[d] >= x[d].first() - 1e-12)
            if (b(-units[d]).castable<parallelBoundary>())
                return -unitXYZ;

    // Binary search. Use <= operator so that at faces points belong to the
    // upper cell. Assign points that are on the upper boundary to the nearest
    // internal cell.

    const labelVector N(this->operator[](0).N());
    const labelVector R
    (
        cmptDivide(N, this->operator[](l).N())
    );

    labelVector ijk;

    for (label d = 0; d < 3; d++)
        ijk[d] =
            Foam::min
            (
                Foam::max
                (
                    findLower(x[d], q[d], 0, lessEqOp<scalar>()),
                    0
                ),
                N[d]-1
            );

    return cmptDivide(ijk, R);
}

}

}
