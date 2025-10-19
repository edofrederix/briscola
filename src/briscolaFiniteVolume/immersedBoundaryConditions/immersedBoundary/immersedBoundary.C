#include "immersedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
void immersedBoundary<MeshType>::setMasks()
{
    // Initialize masks to zero
    mask_ = Zero;
    ghostMask_ = Zero;
    wallAdjMask_ = Zero;

    // Cell centers
    const meshField<vector,MeshType>& cc =
        fvMshMetrics_.cellCenters();

    // Set IB mask fields
    forAllCells(mask_, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);

        if (this->isInside(cc(l,d,ijk)))
        {
            mask_(l,d,ijk) = 1;

            for (int f = 0; f < 6; f++)
            {
                const labelVector fo = faceOffsets[f];

                if (!this->isInside(cc[l][d](ijk+fo)))
                    ghostMask_(l,d,ijk) = 1;
            }
        }
        else
        {
            for (int f = 0; f < 6; f++)
            {
                const labelVector fo = faceOffsets[f];

                if (this->isInside(cc[l][d](ijk+fo)))
                    wallAdjMask_(l,d,i,j,k) = 1;
            }
        }
    }

    // Ghost cells

    forAll(mask_, l)
    {
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
                    if (this->isInside(cc(l,d,ijk+bo)))
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
    // Initialize wall distances to zero
    wallDistAdj_ = Zero;
    wallDistGhost_ = Zero;

    // Cell centers
    const meshField<vector,MeshType>& cc =
        fvMshMetrics_.cellCenters();

    // To-do: write as face loop

    // Set wall distance fields
    forAllCells(mask_, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const faceLabel I = this->fvMsh_.template I<MeshType>(l,d);

        if (!this->isInside(cc(l,d,ijk)))
        {
            const vector c(cc(l,d,ijk));

            for (int f = 0; f < 6; f++)
            {
                const labelVector fo = faceOffsets[f];
                const label fd = f/2;

                if (this->isInside(cc[l][d](ijk+fo)))
                {
                    // Neighbor cell in the IB
                    const vector nb(cc[l][d](ijk+fo));
                    const scalar wd = this->wallDistance(c,nb);
                    const scalar xi = (Foam::mag(c-nb)-wd)/Foam::mag(c-nb);

                    wallDistAdj_[fd](l,d,upperFaceNeighbor(ijk,f)) = xi;
                }
            }
        }
        else
        {
            const vector gc(cc(l,d,ijk));

            for (int f = 0; f < 6; f++)
            {
                const labelVector fo = faceOffsets[f];
                const label fd = f/2;

                if (!this->isInside(cc[l][d](ijk+fo)))
                {
                    // Wall-adjacent cell
                    const vector wa(cc[l][d](ijk+fo));
                    const scalar wd = this->wallDistance(wa,gc);
                    const scalar xi = (Foam::mag(gc-wa)-wd)/Foam::mag(gc-wa);

                    const labelVector upp(upperFaceNeighbor(ijk,f));

                    wallDistGhost_[fd](l,d,upp) = xi;

                    // Second neighbor
                    vector sn(wa + (wa - cc[l][d](ijk)));

                    if
                    (
                           (i + 2*fo.x() >= I.left())
                        && (i + 2*fo.x() <= I.right())
                        && (j + 2*fo.y() >= I.bottom())
                        && (j + 2*fo.y() <= I.top())
                        && (k + 2*fo.z() >= I.aft())
                        && (k + 2*fo.z() <= I.fore())
                    )
                    {
                        sn = cc[l][d](ijk + 2*fo);
                    }

                    const scalar xi2 = Foam::mag(gc-sn)/Foam::mag(gc-wa);

                    neighborDist_[fd](l,d,upp) = xi2;

                    break;
                }
            }
        }
    }

    forAll(wallDistAdj_, fd)
        wallDistAdj_[fd].correctCommsBoundaryConditions();

    forAll(wallDistGhost_, fd)
        wallDistGhost_[fd].correctCommsBoundaryConditions();

    forAll(neighborDist_, fd)
        neighborDist_[fd].correctCommsBoundaryConditions();
}

template<class MeshType>
void immersedBoundary<MeshType>::setMirrorPoints()
{
    // Cell centers
    const meshField<vector,MeshType>& cc =
        fvMshMetrics_.cellCenters();

    // Mesh
    const mesh& msh = fvMsh_.msh();

    scalar tol = 1e-10;

    forAllCells(mirrorPoints_, l, d, i, j, k)
    {
        const labelVector ijk(i,j,k);

        if (ghostMask_(l,d,ijk))
        {
            vector mp = mirrorPoint(cc(l,d,ijk));

            // To-do: what if the domain is not rectangular or not aligned with
            // the coordinate system?

            // Fix situations where the mirror point is just outside of the mesh
            // bounding box due to rounding errors

            for (int fd = 0; fd < 3; fd++)
            {
                if
                (
                    (mp[fd] <= msh[l].boundingBox()[fd*2] + tol)
                 && (mp[fd] >= msh[l].boundingBox()[fd*2] - tol)
                )
                {
                    mp[fd] += tol;
                }

                if
                (
                    (mp[fd] >= msh[l].boundingBox()[fd*2+1] - tol)
                 && (mp[fd] <= msh[l].boundingBox()[fd*2+1] + tol)
                )
                {
                    mp[fd] -= tol;
                }
            }

            mirrorPoints_(l,d,ijk) = mp;
        }
        else
        {
            // The vector (0,0,0) is still a valid coordinate so the
            // mirrorPoints_ field should only be evaluated at ghost cells

            mirrorPoints_(l,d,ijk) = Zero;
        }
    }

    mirrorPoints_.correctCommsBoundaryConditions();
}

// Constructor

template<class MeshType>
immersedBoundary<MeshType>::immersedBoundary
(
    const fvMeshMetrics<MeshType>& metrics,
    const dictionary& dict
)
:
    fvMsh_(metrics.fvMsh()),
    fvMshMetrics_(metrics),
    name_(dict.dictName()),
    shapeOverlap_(false),
    mask_
    (
        "mask",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
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
        true
    ),
    wallAdjMask_
    (
        "ghostMask",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
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
        true
    ),
    wallDistGhost_
    (
        "wallDistGhost",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    neighborDist_
    (
        "neighborDist",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
        true
    ),
    mirrorPoints_
    (
        "mirrorPoints",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true,
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
                immersedBoundaryShape::New
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
    checkOverlap();
}

template<class MeshType>
immersedBoundary<MeshType>::~immersedBoundary()
{}

template<class MeshType>
bool immersedBoundary<MeshType>::isInside(vector xyz) const
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
scalar immersedBoundary<MeshType>::wallDistance(vector c, vector nb) const
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
            << "Neighbor point should be inside the immersed boundary."
            << endl;
        FatalError.exit();
    }

    scalar dist = -1;

    for (int s = 0; s < shapes_.size(); s++)
    {
        if (Foam::mag(dist+1) < 0.01)
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
        || (dist > mag(c-nb)*1.01)
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
) const
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
        if (Foam::mag(dist + 1.0) < 0.01)
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
) const
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
        if (Foam::mag(dist+1) < 0.01)
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

template<class MeshType>
void immersedBoundary<MeshType>::checkOverlap()
{
    // Cell centers
    const meshField<vector,MeshType>& cc =
        fvMshMetrics_.cellCenters();

    forAllCells(cc,l,d,i,j,k)
    {
        Switch inside = false;

        // Check if xyz is inside any of the IB shapes
        for (int s = 0; s < shapes_.size(); s++)
        {
            if(shapes_[s].isInside(cc(l,d,i,j,k)))
            {
                if(inside == true)
                {
                    shapeOverlap_ = true;
                }
                else
                {
                    inside = true;
                }
            }
        }
    }
}

defineTemplateTypeNameAndDebug(colocatedImmersedBoundary, 0);
defineTemplateTypeNameAndDebug(staggeredImmersedBoundary, 0);

template class immersedBoundary<colocated>;
template class immersedBoundary<staggered>;

}

}

}
