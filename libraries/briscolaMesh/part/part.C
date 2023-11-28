#include "part.H"
#include "mesh.H"
#include "SubList.H"

namespace Foam
{

namespace briscola
{

void part::calcPoints(const mesh& msh, const part* l)
{
    points_.clear();
    points_.setSizeFromCells(N_);

    if (l_ == 0)
    {
        const labelVector start = msh.decomp().myBrickPartStart();

        const brick& b = msh.bricks()[msh.decomp().myBrickNum()];
        const labelVector Nb = b.N();

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
            points_(i,j,k) = l->points()(i*R_.x(), j*R_.y(), k*R_.z());
        }
    }
}

void part::calcGhostPoints(const mesh& msh)
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

    // Next, at parallel and periodic boundaries, set the ghost points equal to
    // the neighboring points

    List<const boundary*> boundaryPtrs(0);

    forAll(msh.boundaries(), bi)
    {
        const word type = msh.boundaries()[bi].type();

        if (type == "parallel" || type == "periodic")
            boundaryPtrs.append(&msh.boundaries()[bi]);
    }

    const label Np = boundaryPtrs.size();

    labelList neighborProcNums(Np);

    PtrList<vectorBlock> sendBuffers(Np);
    PtrList<vectorBlock> recvBuffers(Np);

    forAll(boundaryPtrs, bi)
    {
        const boundary& b = *boundaryPtrs[bi];

        const labelVector bo(b.offset());
        const labelTensor T(b.T());

        const labelVector S(pS(bo));
        const labelVector E(pE(bo));
        const labelVector N(pN(bo));

        neighborProcNums[bi] =
            readLabel(b.dict().lookup("neighborProcNum"));

        sendBuffers.set
        (
            bi,
            new vectorBlock(N, Zero)
        );

        recvBuffers.set
        (
            bi,
            new vectorBlock(cmptMag(T.T() & N))
        );

        // For parallel boundaries, send the point coordinate. For periodic
        // boundaries, send the point-to-point distance vector.

        labelVector ijk;

        if (b.type() == "parallel")
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                sendBuffers[bi](ijk-S) = points_(ijk-bo);
            }
        }
        else
        {
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                sendBuffers[bi](ijk-S) = points_(ijk-bo)-points_(ijk);
            }
        }
    }

    // Send/receive

    forAll(boundaryPtrs, bi)
    {
        const label neighborProcNum = neighborProcNums[bi];

        const label tag =
            readLabel(boundaryPtrs[bi]->dict().lookup("tag"));

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
        const boundary& b = *boundaryPtrs[bi];

        const labelVector bo(b.offset());
        const labelTensor T(b.T());

        const labelVector S(pS(bo));
        const labelVector E(pE(bo));

        vectorBlock& recv = recvBuffers[bi];

        recv.transform(T);

        // Set points

        labelVector ijk;

        if (b.type() == "parallel")
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

part::part(const mesh& msh, const part* l)
:
    l_(l == nullptr ? 0 : l->partNum()+1)
{
    if (l == nullptr)
    {
        N_ = msh.decomp().myPartN();
        R_ = zeroXYZ;
    }
    else
    {
        const labelVector P(l->N());

        N_ =
            labelVector
            (
                P.x() < 4 ? P.x() : P.x()/2,
                P.y() < 4 ? P.y() : P.y()/2,
                P.z() < 4 ? P.z() : P.z()/2
            );

        R_ = cmptDivide(P,N_);
    }

    calcPoints(msh, l);
    calcGhostPoints(msh);

    // Check if part directions are rectilinear, which is true if all left,
    // bottom and aft face vectors are aligned

    const scalar tol = 1e-12;

    vector x = zeroXYZ;
    vector y = zeroXYZ;
    vector z = zeroXYZ;

    rectilinear_ = unitXYZ;

    for (int i = 0; i < this->l(); i++)
    for (int j = 0; j < this->m(); j++)
    for (int k = 0; k < this->n(); k++)
    {
        labelVector ijk(i,j,k);

        const vector left =
            0.5
          * (
                (
                    points_(ijk+unitYZ)
                  - points_(ijk)
                )
              ^ (
                    points_(ijk+unitZ)
                  - points_(ijk+unitY)
                )
            );

        const vector bottom =
            0.5
          * (
                (
                    points_(ijk+unitZ)
                  - points_(ijk+unitX)
                )
              ^ (
                    points_(ijk+unitXZ)
                  - points_(ijk)
                )
            );

        const vector aft =
            0.5
          * (
                (
                    points_(ijk+unitXY)
                  - points_(ijk)
                )
              ^ (
                    points_(ijk+unitY)
                  - points_(ijk+unitX)
                )
            );

        if (i == 0 && j == 0 && k == 0)
        {
            x = left;
            y = bottom;
            z = aft;
        }
        else
        {
            if (Foam::mag(x ^ left)/Foam::mag(x) > tol)
            {
                rectilinear_.x() = 0;
            }

            if (Foam::mag(y ^ bottom)/Foam::mag(y) > tol)
            {
                rectilinear_.y() = 0;
            }

            if (Foam::mag(z ^ aft)/Foam::mag(z) > tol)
            {
                rectilinear_.z() = 0;
            }
        }
    }

    // Check if part directions are uniform, which is true if all cell sizes in
    // a direction are the same

    uniform_ = unitXYZ;

    for (int d = 0; d < 3; d++)
    if (rectilinear_[d])
    {
        const vector base = units[d];
        const scalar size = Foam::mag(points_(base) - points_(zeroXYZ));

        for (int i = 1; i < this->N()[d]; i++)
        {
            const scalar next =
                Foam::mag(points_(base*(i+1)) - points_(base*i));

            if (Foam::mag(next-size) > tol)
            {
                uniform_[d] = 0;
                break;
            }
        }
    }
    else
    {
        uniform_[d] = 0;
    }

    // Compute bounding box

    scalarBlock xp(this->N()+unitXYZ);
    scalarBlock yp(this->N()+unitXYZ);
    scalarBlock zp(this->N()+unitXYZ);

    forAllBlock(xp, i, j, k)
    {
        xp(i,j,k) = points_(i,j,k).x();
        yp(i,j,k) = points_(i,j,k).y();
        zp(i,j,k) = points_(i,j,k).z();
    }

    boundingBox_ =
        faceScalar(min(xp), max(xp), min(yp), max(yp), min(zp), max(zp));
}

part::part(const part& l)
:
    l_(l.l_),
    N_(l.N_),
    R_(l.R_),
    points_(l.points_),
    rectilinear_(l.rectilinear_),
    uniform_(l.uniform_),
    boundingBox_(l.boundingBox_)
{}

part::~part()
{}

}

}
