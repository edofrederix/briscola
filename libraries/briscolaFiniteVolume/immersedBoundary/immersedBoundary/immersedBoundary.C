#include "immersedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

using Foam::max;
using Foam::min;
using Foam::sqr;
using Foam::mag;

template<class MeshType>
void immersedBoundary<MeshType>::setMasks()
{
    // Cell centers
    const meshField<vector,MeshType>& CC =
        fvMsh_.metrics<MeshType>().cellCenters();

    // Set IB mask fields
    forAllCells(mask_,l,d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        mask_(l,d,i,j,k) = Zero;
        ghostMask_(l,d,i,j,k) = Zero;
        wallAdjMask_(l,d,i,j,k) = Zero;

        if (this->isInside(CC(l,d,i,j,k)))
        {
            mask_(l,d,i,j,k) = 1;

            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];

                if (!this->isInside(CC[l][d](ijk+fo)))
                {
                    ghostMask_(l,d,i,j,k) = 1;
                }
            }
        }
        else
        {
            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];

                if (this->isInside(CC[l][d](ijk+fo)))
                {
                    wallAdjMask_(l,d,i,j,k) = 1;
                }
            }
        }
    }

    // Ghost cells

    forAll(mask_, l)
    {
        const labelVector N(fvMsh_[l].N());

        forAll(mask_[l], d)
        {
            labelVector bo;

            for (bo.x() = -1; bo.x() <= 1; bo.x()++)
            for (bo.y() = -1; bo.y() <= 1; bo.y()++)
            for (bo.z() = -1; bo.z() <= 1; bo.z()++)
            if (cmptSum(cmptMag(bo)) > 0)
            {
                const labelVector S(mask_.fvMsh().template S<MeshType>(l,d,bo));
                const labelVector E(mask_.fvMsh().template E<MeshType>(l,d,bo));

                labelVector ijk;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    if (this->isInside(CC(l,d,ijk+bo)))
                    {
                        mask_(l,d,ijk+bo) = 1;
                    }
                }
            }
        }
    }

    mask_.correctCommsBoundaryConditions();
    ghostMask_.correctCommsBoundaryConditions();
    wallAdjMask_.correctCommsBoundaryConditions();
}

template<class MeshType>
void immersedBoundary<MeshType>::calculateWallDistances()
{
    // Cell centers
    const meshField<vector,MeshType>& CC =
        fvMsh_.metrics<MeshType>().cellCenters();

    // Set wall distance fields
    forAllCells(mask_,l,d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        wallDistAdj_(l,d,i,j,k) = Zero;
        wallDistGhost_(l,d,i,j,k) = Zero;

        if (!this->isInside(CC(l,d,i,j,k)))
        {
            const vector c(CC(l,d,i,j,k));

            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];

                if (this->isInside(CC[l][d](ijk+fo)))
                {
                    // Neighbor cell in the IB
                    const vector nb(CC[l][d](ijk+fo));
                    const scalar wd = this->wallDistance(c,nb);
                    const scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                    wallDistAdj_(l,d,i,j,k)[dir] = xi;
                }
            }
        }

        if (this->isInside(CC(l,d,i,j,k)))
        {
            const vector gc(CC(l,d,i,j,k));

            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];

                if (!this->isInside(CC[l][d](ijk+fo)))
                {
                    // Wall-adjacent cell
                    const vector wa(CC[l][d](ijk+fo));
                    const scalar wd = this->wallDistance(wa,gc);
                    const scalar xi = (mag(gc-wa)-wd)/mag(gc-wa);

                    wallDistGhost_(l,d,i,j,k)[dir] = xi;

                    // Second neighbor
                    const vector sn(CC[l][d](ijk+2.0*fo));
                    const scalar xi2 = mag(gc-sn)/mag(gc-wa);

                    neighborDist_(l,d,i,j,k)[dir] = xi2;

                    break;
                }
            }
        }
    }

    wallDistAdj_.correctCommsBoundaryConditions();
    wallDistGhost_.correctCommsBoundaryConditions();
    neighborDist_.correctCommsBoundaryConditions();
}

template<class MeshType>
void immersedBoundary<MeshType>::setMirrorPoints()
{
    // Cell centers
    const meshField<vector,MeshType>& CC =
        fvMsh_.metrics<MeshType>().cellCenters();

    // Mesh
    const mesh& msh = fvMsh_.msh();

    scalar tol = 1e-5;

    forAllCells(mirrorPoints_,l,d,i,j,k)
    {
        if (ghostMask_(l,d,i,j,k))
        {
            vector mp = mirrorPoint(CC(l,d,i,j,k));

            // Fix situations where the mirror point is just outside of the mesh
            // bounding box due to rounding errors
            if
            (
                (mp.x() <= msh[l].boundingBox().left() + tol)
                && (mp.x() >= msh[l].boundingBox().left() - tol)
            )
            {
                mp.x() += tol;
            }
            if
            (
                (mp.x() >= msh[l].boundingBox().right() - tol)
                && (mp.x() <= msh[l].boundingBox().right() + tol)
            )
            {
                mp.x() -= tol;
            }
            if
            (
                (mp.y() <= msh[l].boundingBox().bottom() + tol)
                && (mp.y() >= msh[l].boundingBox().bottom() - tol)
            )
            {
                mp.y() += tol;
            }
            if
            (
                (mp.y() >= msh[l].boundingBox().top() - tol)
                && (mp.y() <= msh[l].boundingBox().top() + tol)
            )
            {
                mp.y() -= tol;
            }
            if
            (
                (mp.z() <= msh[l].boundingBox().aft() + tol)
                && (mp.z() >= msh[l].boundingBox().aft() - tol)
            )
            {
                mp.z() += tol;
            }
            if
            (
                (mp.z() >= msh[l].boundingBox().fore() - tol)
                && (mp.z() <= msh[l].boundingBox().fore() + tol)
            )
            {
                mp.z() -= tol;
            }

            mirrorPoints_(l,d,i,j,k) = mp;
        }
        else
        {
            // The vector (0,0,0) is still a valid coordinate
            // so the mirrorPoints_ field should only be evaluated
            // at ghost cells
            mirrorPoints_(l,d,i,j,k) = Zero;
        }
    }

    mirrorPoints_.correctCommsBoundaryConditions();
}

// Constructor

template<class MeshType>
immersedBoundary<MeshType>::immersedBoundary
(
    const fvMesh& fvMsh,
    const dictionary& dict
)
:
    fvMsh_(fvMsh),
    name_(dict.dictName()),
    mask_
    (
        "mask",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        true
    ),
    ghostMask_
    (
        "ghostMask",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        true
    ),
    wallAdjMask_
    (
        "ghostMask",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        true
    ),
    wallDistAdj_
    (
        "wallDistAdj",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        true
    ),
    wallDistGhost_
    (
        "wallDistGhost",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true,
        true
    ),
    neighborDist_
    (
        "neighborDist",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        true
    ),
    mirrorPoints_
    (
        "mirrorPoints",
        fvMsh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        false,
        true
    )
{
    int nShapes = dict.size();

    // Add shapes to IB according to dictionary entries
    for (int e = 0; e < nShapes; e++)
    {
        word entry = dict.toc()[e];

        if (dict.isDict(entry))
        {
            const dictionary& entryDict = dict.subDict(entry);

            shapes_.append
            (
                shape::New
                (
                    entryDict,
                    bool(entryDict.lookupOrDefault("inverted", false))
                )
            );
        }
    }

    setMasks();
    calculateWallDistances();
    setMirrorPoints();
}

template<class MeshType>
immersedBoundary<MeshType>::~immersedBoundary()
{}

template<class MeshType>
bool immersedBoundary<MeshType>::isInside(vector xyz)
{
    // Check if xyz is inside any of the IB shapes
    for (int s = 0; s < shapes_.size(); s++)
    {
        if(shapes_[s].isInside(xyz))
        {
            return true;
        }
    }

    return false;
}

template<class MeshType>
scalar immersedBoundary<MeshType>::wallDistance(vector c, vector nb)
{
    if (this->isInside(c))
    {
        FatalError
            << "Central point should be fluid node."
            << endl;
        FatalError.exit();
    }

    if (!this->isInside(nb))
    {
        FatalError
            << "Neighbor point should be inside the immersed boundary." << c << nb
            << endl;
        FatalError.exit();
    }

    scalar dist = -1;

    for (int s = 0; s < shapes_.size(); s++)
    {
        if (mag(dist+1) < 0.01)
        {
            dist = shapes_[s].wallDistance(c, nb);
        }

        if
        (
               (shapes_[s].wallDistance(c, nb) >= 0)
            && (shapes_[s].wallDistance(c, nb) <= dist)
        )
        {
            dist = shapes_[s].wallDistance(c, nb);
        }
    }

    if
    (
           (dist < 0)
        || (dist > mag(c-nb))
    )
    {
        FatalError
            << "No immersed boundary intersection found."
            << " Distance = " << dist << ", mag(c-nb) = " << mag(c-nb)
            << endl;
        FatalError.exit();
    }

    return dist;
}

template<class MeshType>
scalar immersedBoundary<MeshType>::wallNormalDistance
(
    vector gc
)
{
    if (!this->isInside(gc))
    {
        FatalError
            << "Ghost cell should be inside the immersed boundary."
            << endl;
        FatalError.exit();
    }

    scalar dist = -1;

    for (int s = 0; s < shapes_.size(); s++)
    {
        if (mag(dist+1) < 0.01)
        {
            dist = shapes_[s].wallNormalDistance(gc);
        }

        if
        (
               (shapes_[s].wallNormalDistance(gc) >= 0)
            && (shapes_[s].wallNormalDistance(gc) <= dist)
        )
        {
            dist = shapes_[s].wallNormalDistance(gc);
        }
    }

    if (dist < 0)
    {
        FatalError
            << "No immersed boundary intersection found."
            << " Distance = " << dist << ", for ghost cell = "
            << gc
            << endl;
        FatalError.exit();
    }

    return dist;
}

template<class MeshType>
vector immersedBoundary<MeshType>::mirrorPoint
(
    vector gc
)
{
    if (!this->isInside(gc))
    {
        FatalError
            << "Ghost cell should be inside the immersed boundary."
            << endl;
        FatalError.exit();
    }

    vector mirror = gc;
    scalar dist = -1;

    for (int s = 0; s < shapes_.size(); s++)
    {
        if (mag(dist+1) < 0.01)
        {
            dist = shapes_[s].wallNormalDistance(gc);
            mirror = shapes_[s].mirrorPoint(gc);
        }

        if
        (
               (shapes_[s].wallNormalDistance(gc) >= 0)
            && (shapes_[s].wallNormalDistance(gc) <= dist)
        )
        {
            dist = shapes_[s].wallNormalDistance(gc);
            mirror = shapes_[s].mirrorPoint(gc);
        }
    }

    if (dist < 0)
    {
        FatalError
            << "No immersed boundary intersection found."
            << " Distance = " << dist << ", for ghost cell = "
            << gc
            << endl;
        FatalError.exit();
    }

    return mirror;
}

defineTemplateTypeNameAndDebug(colocatedImmersedBoundary, 0);
defineTemplateTypeNameAndDebug(staggeredImmersedBoundary, 0);

template class immersedBoundary<colocated>;
template class immersedBoundary<staggered>;

}

}

}
