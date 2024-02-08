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

// Constructor

template<class MeshType>
immersedBoundary<MeshType>::immersedBoundary
(
    const fvMesh& fvMsh
)
:
    fvMsh_(fvMsh),
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
    )
{
    const dictionary& IBMdict = fvMsh.meshDict().subDict("ImmersedBoundary");

    int nEntries = IBMdict.size();

    // Add shapes to IB according to dictionary entries
    for (int e = 0; e < nEntries; e++)
    {
        word entry = IBMdict.toc()[e];

        if (IBMdict.isDict(entry))
        {
            const dictionary& entryDict = IBMdict.subDict(entry);

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

    // Cell centers
    const meshField<vector,MeshType>& CC =
        fvMsh_.metrics<MeshType>().cellCenters();

    // Set IB mask fields
    forAllCells(mask_,l,d,i,j,k)
    {
        const labelVector ijk(i,j,k);

        mask_(l,d,i,j,k) = 0;
        ghostMask_(l,d,i,j,k) = 0;
        wallAdjMask_(l,d,i,j,k) = 0;

        if (this->isInside(CC(l,d,i,j,k)))
        {
            mask_(l,d,i,j,k) = 1;

            for (int dir = 0; dir < 6; dir++)
            {
                const labelVector fo = faceOffsets[dir];

                if
                (
                    !this->isInside(CC[l][d](ijk+fo))
                )
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

                if
                (
                    this->isInside(CC[l][d](ijk+fo))
                )
                {
                    wallAdjMask_(l,d,i,j,k) = 1;
                }
            }
        }
    }

    mask_.correctParallelBoundaryConditions();
    ghostMask_.correctParallelBoundaryConditions();
    wallAdjMask_.correctParallelBoundaryConditions();
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
            << "Neighbor point should be inside the immersed boundary."
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
