#include "rectilinearMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "parallelBoundary.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(rectilinearMesh, 0);
addToRunTimeSelectionTable(mesh, rectilinearMesh, dictionary);

void rectilinearMesh::setMetrics(const label l)
{
    const level& lvl = this->operator[](l);
    const levelPoints& points = lvl.points();

    // Base vectors

    if (l == 0)
    {
        vectorList base(3);

        for (int d = 0; d < 3; d++)
        {
            base[d] = points(units[d]) - points(zeroXYZ);
            base[d] /= Foam::mag(base[d]);
        }

        base_ = tensor(base[0], base[1], base[2]);

        if
        (
            base_.x() != base[0]
         || base_.y() != base[1]
         || base_.z() != base[2]
        )
        {
            FatalErrorInFunction
                << "Base tensor does not match base vectors" << endl
                << abort(FatalError);
        }
    }

    // Local cell sizes and points

    localCellSizesData_[l].clear();
    localPointsData_[l].clear();
    localCellSizesData_[l].setSize(3);
    localPointsData_[l].setSize(3);

    for (int d = 0; d < 3; d++)
    {
        if (lvl.empty())
        {
            localCellSizesData_[l].set(d, new scalarList(0));
            localPointsData_[l].set(d, new scalarList(0));
        }
        else
        {
            localCellSizesData_[l].set(d, new scalarList(lvl.N()[d]+2));
            localPointsData_[l].set(d, new scalarList(lvl.N()[d]+3));
        }

        scalarList& localCellSizesData = localCellSizesData_[l][d];
        scalarList& localPointsData = localPointsData_[l][d];

        const labelVector dir = units[d];

        const vector base2 =
            d == 0 ? base_.x() : d == 1 ? base_.y() : base_.z();

        forAll(localPointsData, i)
        {
            localPointsData[i] = trimPrecision(points(dir*(i-1)) & base2);
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

    localCellSizes_[l].clear();
    localPoints_[l].clear();
    localCellSizes_[l].setSize(3);
    localPoints_[l].setSize(3);

    for (int d = 0; d < 3; d++)
    {
        if (lvl.empty())
        {
            localCellSizes_[l].set
            (
                d,
                new PartialList<scalar>(localCellSizesData_[l][d], 0)
            );

            localPoints_[l].set
            (
                d,
                new PartialList<scalar>(localPointsData_[l][d], 0)
            );
        }
        else
        {
            localCellSizes_[l].set
            (
                d,
                new PartialList<scalar>
                (
                    localCellSizesData_[l][d],
                    localCellSizesData_[l][d].size()-2,
                    1
                )
            );

            localPoints_[l].set
            (
                d,
                new PartialList<scalar>
                (
                    localPointsData_[l][d],
                    localPointsData_[l][d].size()-2,
                    1
                )
            );
        }
    }

    // Global cell sizes

    const decompositionMap& map = lvl.decomp().map();

    globalCellSizesData_[l].clear();
    globalPointsData_[l].clear();
    globalCellSizesData_[l].setSize(3);
    globalPointsData_[l].setSize(3);

    PtrList<labelList> starts(3);

    for (int d = 0; d < 3; d++)
    {
        if (lvl.empty())
        {
            globalCellSizesData_[l].set(d, new scalarList(0));
            globalPointsData_[l].set(d, new scalarList(0));
        }
        else
        {
            globalCellSizesData_[l].set(d, new scalarList(this->N()[d]+2));
            globalPointsData_[l].set(d, new scalarList(this->N()[d]+3));
        }

        scalarList& globalCellSizes = globalCellSizesData_[l][d];
        scalarList& globalPoints = globalPointsData_[l][d];

        const labelVector base2 = units[d];

        label cursor = 0;

        starts.set(d, new labelList(map.shape()[d]));

        for (int i = 0; i < map.shape()[d]; i++)
        {
            const label proc = map(base2*i);

            if (proc == Pstream::myProcNo())
            {
                const scalarList& localCellSizesData =
                    localCellSizesData_[l][d];

                const scalarList& localPointsData =
                    localPointsData_[l][d];

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
                        Pstream::masterNo(),
                        0,
                        UPstream::msgType(),
                        lvl.comms()
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
                    proc,
                    0,
                    UPstream::msgType(),
                    lvl.comms()
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

        Pstream::scatter(globalCellSizes, Pstream::msgType(), lvl.comms());
        Pstream::scatter(globalPoints, Pstream::msgType(), lvl.comms());

        if (Pstream::master())
        {
            for (int proc = 0; proc < Pstream::nProcs(); proc++)
            if (proc != Pstream::myProcNo())
            {
                OPstream send
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    0,
                    UPstream::msgType(),
                    lvl.comms()
                );

                send << starts[d];
            }
        }
        else
        {
            IPstream recv
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo(),
                0,
                UPstream::msgType(),
                lvl.comms()
            );

            recv >> starts[d];
        }
    }

    // Create global lists without ghosts as sub-lists

    globalCellSizes_[l].clear();
    globalPoints_[l].clear();
    globalCellSizes_[l].setSize(3);
    globalPoints_[l].setSize(3);

    for (int d = 0; d < 3; d++)
    {
        if (lvl.empty())
        {
            globalCellSizes_[l].set
            (
                d,
                new PartialList<scalar>(globalCellSizesData_[l][d], 0)
            );

            globalPoints_[l].set
            (
                d,
                new PartialList<scalar>(globalPointsData_[l][d], 0)
            );
        }
        else
        {
            globalCellSizes_[l].set
            (
                d,
                new PartialList<scalar>
                (
                    globalCellSizesData_[l][d],
                    globalCellSizesData_[l][d].size()-2,
                    1
                )
            );

            globalPoints_[l].set
            (
                d,
                new PartialList<scalar>
                (
                    globalPointsData_[l][d],
                    globalPointsData_[l][d].size()-2,
                    1
                )
            );
        }
    }

    // Set global starts

    globalStarts_[l].setSize(Pstream::nProcs());

    forAllBlock(map, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const label proc = map(ijk);

        globalStarts_[l][proc] = Zero;

        for (int d = 0; d < 3; d++)
            globalStarts_[l][proc][d] = starts[d][ijk[d]];
    }

    // Check if global mesh directions are uniform, only on the first level

    if (l == 0)
    {
        globalUniform_ = labelVector(1,1,1);

        const scalar tol = 1e-7;

        forAll(globalCellSizes_[l], d)
        {
            scalar maxCellSize = max(globalCellSizes_[l][d]);
            scalar minCellSize = min(globalCellSizes_[l][d]);

            if (maxCellSize - minCellSize > tol)
            {
                globalUniform_[d] = 0;
            }
        }
    }
}

void rectilinearMesh::setMetrics()
{
    const label nLevels = this->size();

    localCellSizesData_.setSize(nLevels);
    localPointsData_.setSize(nLevels);
    globalCellSizesData_.setSize(nLevels);
    globalPointsData_.setSize(nLevels);
    localCellSizes_.setSize(nLevels);
    localPoints_.setSize(nLevels);
    globalCellSizes_.setSize(nLevels);
    globalPoints_.setSize(nLevels);

    globalStarts_.setSize(nLevels);

    forAll(*this, l)
    {
        localCellSizesData_.set(l, new FastPtrList<scalarList>());
        localPointsData_.set(l, new FastPtrList<scalarList>());
        globalCellSizesData_.set(l, new FastPtrList<scalarList>());
        globalPointsData_.set(l, new FastPtrList<scalarList>());

        localCellSizes_.set(l, new FastPtrList<PartialList<scalar>>());
        localPoints_.set(l, new FastPtrList<PartialList<scalar>>());
        globalCellSizes_.set(l, new FastPtrList<PartialList<scalar>>());
        globalPoints_.set(l, new FastPtrList<PartialList<scalar>>());

        setMetrics(l);
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
    base_(msh.base_)
{
    setMetrics();
}

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

    const FastPtrList<PartialList<scalar>>& x = localPoints_[l];

    const level& lvl = this->operator[](l);

    // On empty levels there's nothing to search

    if (lvl.empty())
        return -unitXYZ;

    // Check if the point is outside the domain boundaries

    for (label d = 0; d < 3; d++)
        if (q[d] < x[d].first() || q[d] > x[d].last() + 1e-12)
            return -unitXYZ;

    // Check if the point is just above/below the upper/lower boundary. If so,
    // the point is included if the boundary is not be a parallel one. If the
    // point is exactly on the lower boundary it must be included anyway.

    for (label d = 0; d < 3; d++)
        if (q[d] >= x[d].last() && q[d] <= x[d].last() + 1e-12)
            if (lvl.boundaries().find(units[d]).castable<parallelBoundary>())
                return -unitXYZ;

    for (label d = 0; d < 3; d++)
        if (q[d] < x[d].first() && q[d] >= x[d].first() - 1e-12)
            if (lvl.boundaries().find(-units[d]).castable<parallelBoundary>())
                return -unitXYZ;

    // Binary search. Use <= operator so that at faces points belong to the
    // upper cell. Assign points that are on the upper boundary to the nearest
    // internal cell.

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
                lvl.N()[d]-1
            );

    return ijk;
}

}

}
