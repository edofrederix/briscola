#include "partPoints.H"
#include "mesh.H"

#include "parallelBoundary.H"
#include "periodicBoundary.H"

namespace Foam
{

namespace briscola
{

partPoints::partPoints(const mesh& msh)
:
    vectorBlock(),
    msh_(msh)
{}

partPoints::partPoints(const partPoints& points)
:
    vectorBlock(points),
    msh_(points.msh_)
{}

partPoints::~partPoints()
{}

void partPoints::clear()
{
    vectorBlock::clear();
}

void partPoints::setSizeFromCells(const labelVector& size)
{
    // On each face of the block a layer of ghost points, plus one additional
    // layer because these are points, not cells.

    vectorBlock& points = *this;

    points.setSize(size + 3*unitXYZ);
    points = Zero;
}

void partPoints::calcGhostPoints()
{
    vectorBlock& points = *this;

    const labelVector N = points.shape();

    // Project boundary vertices along vertex lines

    for (int fi = 0; fi < 6; fi++)
    {
        const labelVector bo(faceOffsets[fi]);

        const labelVector S(start(bo,N));
        const labelVector E(end(bo,N));

        labelVector ijk;

        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            points(ijk + bo) = 2.0*points(ijk) - points(ijk-bo);
        }
    }

    // Also set edge and vertex vertices on structured meshes only

    if (msh_.structured())
    {
        for (int ei = 0; ei < 12; ei++)
        {
            const labelVector bo(edgeOffsets[ei]);
            const label dir = ei/4;

            const labelVector S(start(bo,N));
            const labelVector E(end(bo,N));

            labelVector ijk;

            const labelVector bo1(dir == 0 ? bo.y()*unitY : bo.x()*unitX);
            const labelVector bo2(dir == 2 ? bo.y()*unitY : bo.z()*unitZ);

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                points(ijk+bo) =
                    points(ijk+bo1) + points(ijk+bo2) - points(ijk);
            }
        }

        for (int vi = 0; vi < 8; vi++)
        {
            const labelVector bo(vertexOffsets[vi]);

            const labelVector S(start(bo,N));
            const labelVector E(end(bo,N));

            labelVector ijk;

            const labelVector bo1(bo.x()*unitX);
            const labelVector bo2(bo.y()*unitY);
            const labelVector bo3(bo.z()*unitZ);

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                points(ijk+bo) =
                    points(ijk+bo1)
                  + points(ijk+bo2)
                  + points(ijk+bo3)
                  - points(ijk)*2.0;
            }
        }
    }

    // Next, at parallel and periodic boundaries, set the ghost points equal to
    // the neighboring points

    List<const boundary*> boundaryPtrs(0);

    for (int deg = 1; deg < 3; deg++)
        forAll(msh_.boundaries(), bi)
            if (msh_.boundaries()[bi].castable<parallelBoundary>())
                if (msh_.boundaries()[bi].offsetDegree() == deg)
                    boundaryPtrs.append(&msh_.boundaries()[bi]);

    const label Np = boundaryPtrs.size();

    PtrList<vectorBlock> sendBuffers(Np);
    PtrList<vectorBlock> recvBuffers(Np);

    forAll(boundaryPtrs, bi)
    {
        const parallelBoundary& b =
            boundaryPtrs[bi]->cast<parallelBoundary>();

        const labelVector bo(b.offset());
        const labelTensor T(b.T());

        labelVector S(start(bo,N));
        labelVector E(end(bo,N));

        if (msh_.structured())
        {
            labelVector expansion(unitXYZ);

            for (int dir = 0; dir < 3; dir++)
                if(bo[dir] != 0)
                    expansion[dir] = 0;

            S -= expansion;
            E += expansion;
        }

        const labelVector shape(E-S);

        sendBuffers.set
        (
            bi,
            new vectorBlock(shape, Zero)
        );

        recvBuffers.set
        (
            bi,
            new vectorBlock(cmptMag(T.T() & shape))
        );

        // For parallel boundaries, send the point coordinate. For periodic
        // boundaries, send the point-to-point distance vector.

        labelVector ijk;

        if (!b.castable<periodicBoundary>())
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                sendBuffers[bi](ijk-S) = points(ijk-bo);
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                sendBuffers[bi](ijk-S) = points(ijk-bo)-points(ijk);
            }
        }
    }

    // Send/receive

    forAll(boundaryPtrs, bi)
    {
        const parallelBoundary& b =
            boundaryPtrs[bi]->cast<parallelBoundary>();

        const label neighborProcNum = b.neighborProcNum();
        const label tag = b.tag();

        // Send/receive

        if (neighborProcNum != Pstream::myProcNo())
        {
            UIPstream::read
            (
                Pstream::commsTypes::nonBlocking,
                neighborProcNum,
                reinterpret_cast<char*>(recvBuffers[bi].begin()),
                recvBuffers[bi].byteSize(),
                tag,
                UPstream::worldComm
            );

            UOPstream::write
            (
                Pstream::commsTypes::nonBlocking,
                neighborProcNum,
                reinterpret_cast<char*>(sendBuffers[bi].begin()),
                sendBuffers[bi].byteSize(),
                tag,
                UPstream::worldComm
            );
        }
        else
        {
            // On the same processor the boundary must be periodic, and its
            // neighbor is assumed to be on the opposing face/edge/vertex. Copy
            // buffers directly.

            const boundary& b = *boundaryPtrs[bi];

            const labelVector bo(b.offset());

            forAll(boundaryPtrs, bj)
            if (boundaryPtrs[bj]->offset() == -bo)
            {
                recvBuffers[bi] = sendBuffers[bj];
                break;
            }
        }
    }

    UPstream::waitRequests();

    // Unpack

    forAll(boundaryPtrs, bi)
    {
        const parallelBoundary& b =
            boundaryPtrs[bi]->cast<parallelBoundary>();

        const labelVector bo(b.offset());
        const labelTensor T(b.T());

        labelVector S(start(bo,N));
        labelVector E(end(bo,N));

        if (msh_.structured())
        {
            labelVector expansion(unitXYZ);

            for (int dir = 0; dir < 3; dir++)
                if(bo[dir] != 0)
                    expansion[dir] = 0;

            S -= expansion;
            E += expansion;
        }

        vectorBlock& recv = recvBuffers[bi];

        recv.transform(T);

        // Set points

        labelVector ijk;

        if (!b.castable<periodicBoundary>())
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                points(ijk+bo) = recv(ijk-S);
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                // For periodic points we project the neighbor's point-to-point
                // distance along the current cell's point-to-point vector. This
                // is only correct for orthogonal periodicity.

                const vector n = normalised(points(ijk) - points(ijk-bo));

                points(ijk+bo) = points(ijk) + n*mag(recv(ijk-S));
            }
        }
    }
}

void partPoints::clean()
{
    vectorBlock& points = *this;

    forAllBlock(points, i, j, k)
        points(i,j,k) = briscola::clean(points(i,j,k));
}

}

}
