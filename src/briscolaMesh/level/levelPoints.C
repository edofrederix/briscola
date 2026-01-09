#include "levelPoints.H"
#include "mesh.H"

#include "parallelBoundary.H"
#include "periodicBoundary.H"

#include "hexa.H"

namespace Foam
{

namespace briscola
{

levelPoints::levelPoints(const level& lvl)
:
    vectorBlock(),
    lvl_(lvl)
{
    setSizeFromCells(lvl_.N());

    // Reference needs to be of levelPoints type to account for ghost point
    // padding

    levelPoints& points = *this;

    if (lvl_.levelNum() == 0)
    {
        const brick& b = lvl_.msh().bricks()[lvl_.decomp().myBrickNum()];
        const labelVector Nb = b.N();

        const labelVector start =
            cmptMultiply
            (
                cmptDivide(Nb, lvl_.decomp().myBrickDecomp()),
                lvl_.decomp().myBrickPart()
            );

        // Normalized coordinates of the points in this processor's part of the
        // brick

        scalarField xi(lvl_.N().x()+1);
        scalarField eta(lvl_.N().y()+1);
        scalarField zeta(lvl_.N().z()+1);

        forAll(xi, i)
            xi[i] = scalar(start.x()+i)/Nb.x();

        forAll(eta, j)
            eta[j] = scalar(start.y()+j)/Nb.y();

        forAll(zeta, k)
            zeta[k] = scalar(start.z()+k)/Nb.z();

        const vectorBlock tfi(b.TFI(xi,eta,zeta));

        forAllBlock(tfi, i, j, k)
            points(i,j,k) = tfi(i,j,k);
    }
    else
    {
        // Sample parent level points

        const level& parent = lvl_.parent();
        const labelVector R = lvl_.R();

        forAllBlock(points, i, j, k)
            points(i,j,k) = parent.points()(i*R.x(), j*R.y(), k*R.z());
    }

    // Set ghost points

    calcGhostPoints();

    // Remove round-off errors

    clean();
}

levelPoints::levelPoints(const levelPoints& points)
:
    vectorBlock(points),
    lvl_(points.lvl_)
{}

levelPoints::levelPoints(const levelPoints& points, const level& lvl)
:
    vectorBlock(points),
    lvl_(lvl)
{}

levelPoints::~levelPoints()
{}

void levelPoints::clear()
{
    vectorBlock::clear();
}

void levelPoints::setSizeFromCells(const labelVector& size)
{
    // On each face of the block a layer of ghost points, plus one additional
    // layer because these are points, not cells.

    vectorBlock& points = *this;

    points.setSize(size + 3*unitXYZ);
    points = Zero;
}

void levelPoints::calcGhostPoints()
{
    vectorBlock& points = *this;

    const labelVector N = points.shape();

    // Project boundary vertices along the vector that connects the boundary
    // vertex to its internal vertex neighbor. This may cause invalid ghost cell
    // hexahedrals especially at coarser grid levels. Once an invalid hexahedral
    // is encountered, the process is repeated across all processors with half
    // the projection distance, until a valid mesh is found.

    scalar dist(1.0);
    label nIter(0);

    while (true)
    {
        nIter++;

        // Face ghost points

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
                points(ijk+bo) =
                    points(ijk) + dist*(points(ijk) - points(ijk-bo));
            }
        }

        // Edge ghost points

        for (int ei = 0; ei < 12; ei++)
        {
            const labelVector bo(edgeOffsets[ei]);
            const label dir = ei/4;

            const labelVector S(start(bo,N));
            const labelVector E(end(bo,N));

            const labelVector bo1(dir == 0 ? bo.y()*unitY : bo.x()*unitX);
            const labelVector bo2(dir == 2 ? bo.y()*unitY : bo.z()*unitZ);

            labelVector ijk;
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                // Sum of the two outward ghost vertex vectors at the edge

                points(ijk+bo) =
                    points(ijk+bo1) + points(ijk+bo2) - points(ijk);
            }
        }

        // Vertex ghost point

        for (int vi = 0; vi < 8; vi++)
        {
            const labelVector bo(vertexOffsets[vi]);

            const labelVector ijk(start(bo,N));

            const labelVector bo1(bo.x()*unitX);
            const labelVector bo2(bo.y()*unitY);
            const labelVector bo3(bo.z()*unitZ);

            // Sum of the three outward vertex vectors

            points(ijk+bo) =
                points(ijk+bo1)
              + points(ijk+bo2)
              + points(ijk+bo3)
              - points(ijk)*2.0;
        }

        // Check cell validity

        bool valid = true;

        const labelVector S(start(-unitXYZ,N)-unitXYZ);
        const labelVector E(end(unitXYZ,N));

        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        if (valid)
        {
            const hexa h
            (
                points(ijk),
                points(ijk + unitX),
                points(ijk + unitY),
                points(ijk + unitXY),
                points(ijk + unitZ),
                points(ijk + unitXZ),
                points(ijk + unitYZ),
                points(ijk + unitXYZ)
            );

            if (!h.valid())
                valid = false;
        }

        if (returnReduce(valid, andOp<bool>())) break;

        dist = dist/2.0;

        if (nIter == 1000)
            FatalErrorInFunction
                << "Ghost point calculation failed in " << nIter
                << " iterations" << endl << abort(FatalError);
    }

    // Next, at parallel and periodic boundaries, set the ghost points equal to
    // the neighboring points

    List<const boundary*> boundaryPtrs(0);

    for (int deg = 1; deg < 3; deg++)
        forAll(lvl_.boundaries(), bi)
            if (lvl_.boundaries()[bi].castable<parallelBoundary>())
                if (lvl_.boundaries()[bi].offsetDegree() == deg)
                    boundaryPtrs.append(&lvl_.boundaries()[bi]);

    const label Np = boundaryPtrs.size();

    FastPtrList<vectorBlock> sendBuffers(Np);
    FastPtrList<vectorBlock> recvBuffers(Np);

    forAll(boundaryPtrs, bi)
    {
        const parallelBoundary& b =
            boundaryPtrs[bi]->cast<parallelBoundary>();

        const labelVector bo(b.offset());
        const labelTensor T(b.T());

        labelVector S(start(bo,N));
        labelVector E(end(bo,N));

        if (lvl_.msh().structured())
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

            const boundary& b2 = *boundaryPtrs[bi];

            const labelVector bo(b2.offset());

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

        if (lvl_.msh().structured())
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

void levelPoints::clean()
{
    vectorBlock& points = *this;

    forAllBlock(points, i, j, k)
        points(i,j,k) = briscola::trimPrecision(points(i,j,k));
}

}

}
