#include "partLevel.H"
#include "mesh.H"
#include "SubList.H"

namespace Foam
{

namespace briscola
{

void partLevel::calcPoints()
{
    points_.clear();
    points_.setSizeFromCells(N_);

    if (l_ == 0)
    {
        const labelVector start = msh_.decomp().myBrickPartStart();

        const brick& b = msh_.decomp().myBrick();
        const labelVector Nb = msh_.decomp().myBrickN();

        // Normalized coordinates of the points in this processor's part of the
        // brick

        scalarField xi(N_.x()+1);
        scalarField eta(N_.y()+1);
        scalarField zeta(N_.z()+1);

        forAll(xi, i)
            xi[i] = scalar(start.x()+i)/Nb.x();

        forAll(eta, j)
            eta[j] = scalar(start.y()+j)/Nb.y();

        forAll(zeta, k)
            zeta[k] = scalar(start.z()+k)/Nb.z();

        const vectorBlock tfi(b.TFI(xi,eta,zeta));

        forAllBlock(tfi, i, j, k)
        {
            points_(i,j,k) = tfi(i,j,k);
        }
    }
    else
    {
        forAllBlock(points_, i, j, k)
        {
            points_(i,j,k) = parentPtr_->points()(i*R_.x(), j*R_.y(), k*R_.z());
        }
    }
}

void partLevel::calcGhostPoints()
{
    // First, project all inner points along the point-to-point vector

    const labelVector N = points_.shape();

    labelVector bo;

    for (bo.x() = -1; bo.x() <= 1; bo.x()++)
    for (bo.y() = -1; bo.y() <= 1; bo.y()++)
    for (bo.z() = -1; bo.z() <= 1; bo.z()++)
    if (cmptSum(cmptMag(bo)) > 0)
    {
        const labelVector s
        (
            bo.x() == 0 ? 0 : (bo.x()+1)/2*(N.x()-1),
            bo.y() == 0 ? 0 : (bo.y()+1)/2*(N.y()-1),
            bo.z() == 0 ? 0 : (bo.z()+1)/2*(N.z()-1)
        );

        const labelVector e
        (
            s.x() + (bo.x() == 0 ? N.x()-1 : 0),
            s.y() + (bo.y() == 0 ? N.y()-1 : 0),
            s.z() + (bo.z() == 0 ? N.z()-1 : 0)
        );

        for (label i = s.x(); i <= e.x(); i++)
        for (label j = s.y(); j <= e.y(); j++)
        for (label k = s.z(); k <= e.z(); k++)
        {
            const labelVector ijk(i,j,k);

            points_(ijk+bo) = 2*points_(ijk) - points_(ijk-bo);
        }
    }

    // Next, at parallel and periodic part patches, set the ghost points equal
    // to the neighboring points

    List<const partPatch*> partPatchPtrs(0);

    forAll(msh_.partPatches(), patchi)
    {
        const word type = msh_.partPatches()[patchi].type();

        if (type == "parallel" || type == "periodic")
            partPatchPtrs.append(&msh_.partPatches()[patchi]);
    }

    const label Np = partPatchPtrs.size();

    labelList neighborProcNums(Np);

    PtrList<vectorBlock> sendBuffers(Np);
    PtrList<vectorBlock> recvBuffers(Np);

    forAll(partPatchPtrs, patchi)
    {
        const partPatch& p = *partPatchPtrs[patchi];

        const labelVector bo(p.boundaryOffset());
        const labelTensor T(p.T());

        const labelVector S(pointBoundaryStart(bo));
        const labelVector E(pointBoundaryEnd(bo));

        neighborProcNums[patchi] =
            readLabel(p.dict().lookup("neighborProcNum"));

        sendBuffers.set
        (
            patchi,
            new vectorBlock(N, Zero)
        );

        recvBuffers.set
        (
            patchi,
            new vectorBlock(cmptMag(T.T() & N))
        );

        // For parallel patches, send the point coordinate. For periodic
        // patches, send the point-to-point distance vector.

        labelVector ijk;

        if (p.type() == "parallel")
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                sendBuffers[patchi](ijk-S) = points_(ijk-bo);
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                sendBuffers[patchi](ijk-S) = points_(ijk-bo)-points_(ijk);
            }
        }
    }

    // Send/receive

    forAll(partPatchPtrs, patchi)
    {
        const label neighborProcNum = neighborProcNums[patchi];

        const label tag =
            readLabel(partPatchPtrs[patchi]->dict().lookup("tag"));

        // Send/receive

        if (neighborProcNum != Pstream::myProcNo())
        {
            UIPstream::read
            (
                Pstream::commsTypes::nonBlocking,
                neighborProcNum,
                reinterpret_cast<char*>(recvBuffers[patchi].begin()),
                recvBuffers[patchi].byteSize(),
                tag,
                UPstream::worldComm
            );

            UOPstream::write
            (
                Pstream::commsTypes::nonBlocking,
                neighborProcNum,
                reinterpret_cast<char*>(sendBuffers[patchi].begin()),
                sendBuffers[patchi].byteSize(),
                tag,
                UPstream::worldComm
            );
        }
        else
        {
            // On the same processor the patch part must be periodic, and its
            // neighbor is assumed to be on the opposing face/edge/vertex. Copy
            // buffers directly.

            const partPatch& p = *partPatchPtrs[patchi];

            const labelVector bo(p.boundaryOffset());

            forAll(partPatchPtrs, patchj)
            if (partPatchPtrs[patchj]->boundaryOffset() == -bo)
            {
                recvBuffers[patchi] = sendBuffers[patchj];
                break;
            }
        }
    }

    UPstream::waitRequests();

    // Unpack

    forAll(partPatchPtrs, patchi)
    {
        const partPatch& p = *partPatchPtrs[patchi];

        const labelVector bo(p.boundaryOffset());
        const labelTensor T(p.T());

        const labelVector S(pointBoundaryStart(bo));
        const labelVector E(pointBoundaryEnd(bo));

        vectorBlock& recv = recvBuffers[patchi];

        recv.transform(T);

        // Set points

        labelVector ijk;

        if (p.type() == "parallel")
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                points_(ijk+bo) = recv(ijk-S);
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

                const vector n = normalised(points_(ijk) - points_(ijk-bo));

                points_(ijk+bo) = points_(ijk) + n*mag(recv(ijk-S));
            }
        }
    }
}

partLevel::partLevel(const mesh& msh, const partLevel* l)
:
    msh_(msh),
    parentPtr_(l),
    l_(l == nullptr ? 0 : l->partLevelNum()+1)
{
    if (parentPtr_ == nullptr)
    {
        N_ = msh_.N();
        R_ = zeroXYZ;
    }
    else
    {
        const labelVector P(parentPtr_->N());

        N_ =
            labelVector
            (
                P.x() < 4 ? P.x() : P.x()/2,
                P.y() < 4 ? P.y() : P.y()/2,
                P.z() < 4 ? P.z() : P.z()/2
            );

        R_ = cmptDivide(P,N_);
    }

    calcPoints();
    calcGhostPoints();
}

partLevel::~partLevel()
{}

}

}
