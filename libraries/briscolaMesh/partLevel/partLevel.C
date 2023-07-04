#include "partLevel.H"
#include "mesh.H"
#include "SubList.H"

namespace Foam
{

namespace briscola
{

void partLevel::calcPoints(const mesh& msh, const partLevel* l)
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

void partLevel::calcGhostPoints(const mesh& msh)
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

    forAll(msh.partPatches(), patchi)
    {
        const word type = msh.partPatches()[patchi].type();

        if (type == "parallel" || type == "periodic")
            partPatchPtrs.append(&msh.partPatches()[patchi]);
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
        const labelVector N(pointBoundaryN(bo));

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
    l_(l == nullptr ? 0 : l->partLevelNum()+1)
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

partLevel::partLevel(const partLevel& l)
:
    l_(l.l_),
    N_(l.N_),
    R_(l.R_),
    points_(l.points_),
    rectilinear_(l.rectilinear_),
    uniform_(l.uniform_),
    boundingBox_(l.boundingBox_)
{}

partLevel::~partLevel()
{}

bool partLevel::insideCell
(
    const vector& point,
    const label i,
    const label j,
    const label k
) const
{
    // Bounding box check

    const vectorBlock points(points_.cellPoints(i,j,k));

    const scalarBlock x(points & vector(1,0,0));
    const scalarBlock y(points & vector(0,1,0));
    const scalarBlock z(points & vector(0,0,1));

    const scalar xMin(min(x));
    const scalar xMax(max(x));
    const scalar yMin(min(y));
    const scalar yMax(max(y));
    const scalar zMin(min(z));
    const scalar zMax(max(z));

    const scalar tol = 1e-14;

    if
    (
        point.x() < xMin-tol || point.x() > xMax+tol
     || point.y() < yMin-tol || point.y() > yMax+tol
     || point.z() < zMin-tol || point.z() > zMax+tol
    )
    {
        return false;
    }

    // Point is at least close. Try to find the parameters (u,v,w) in the
    // parametrized function that tri-linearly interpolates the cell. This
    // function properly accounts for faces being doubly ruled surfaces. If
    // (u,v,w) is in the unit cube, the point is in the cell. Find (u,v,w) with
    // Newton's method.

    const vector P(point - points(0));
    const vectorBlock D(points - points(0));

    const vector F1(D(1,0,0));
    const vector F2(D(0,1,0));
    const vector F3(D(0,0,1));

    const vector F4(D(1,1,0)-F1-F2);
    const vector F5(D(1,0,1)-F1-F3);
    const vector F6(D(0,1,1)-F2-F3);

    const vector F7(D(1,1,1)-F1-F2-F3-F4-F5-F6);

    // Initial guess (exact for rectangular cells)

    vector u
    (
        (P & D(1,0,0))/Foam::magSqr(D(1,0,0)),
        (P & D(0,1,0))/Foam::magSqr(D(0,1,0)),
        (P & D(0,0,1))/Foam::magSqr(D(0,0,1))
    );

    const label maxIter = 100;

    vector du(vector::one);
    label iter = 0;

    while (mag(du) > tol && iter < maxIter)
    {
        const tensor dfdu
        (
            F1 + u.y()*F4 + u.z()*F5 + u.y()*u.z()*F7,
            F2 + u.x()*F4 + u.z()*F6 + u.x()*u.z()*F7,
            F3 + u.x()*F5 + u.y()*F6 + u.x()*u.y()*F7
        );

        const vector f
        (
            u.x()*F1
          + u.y()*F2
          + u.z()*F3
          + u.x()*u.y()*F4
          + u.x()*u.z()*F5
          + u.y()*u.z()*F6
          + u.x()*u.y()*u.z()*F7
          - P
        );

        du = inv(dfdu.T()) & f;
        u -= du;

        iter++;
    }

    if (iter == maxIter)
    {
        FatalErrorInFunction
            << "Could not determine if point " << point << " is in cell "
            << labelVector(i,j,k) << endl << abort(FatalError);
    }

    // By convention, a point on a face belongs to the upper cell in the local
    // coordinate direction, unless the point is on an upper boundary.

    return
        u.x() >= 0 && (u.x() < 1 || (u.x() <= 1+tol && i == this->l()-1))
     && u.y() >= 0 && (u.y() < 1 || (u.y() <= 1+tol && j == this->m()-1))
     && u.z() >= 0 && (u.z() < 1 || (u.z() <= 1+tol && k == this->n()-1));
}

}

}
