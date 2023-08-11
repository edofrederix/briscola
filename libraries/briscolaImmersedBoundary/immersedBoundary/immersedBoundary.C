#include "immersedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace ibm
{

using Foam::max;
using Foam::min;

// Constructor

template<class MeshType>
immersedBoundary<MeshType>::immersedBoundary
(
    IOdictionary solverDict,
    const fvMesh& fvMsh
)
:
    fvMsh_(fvMsh),
    mask_("mask", fvMsh_),
    wallAdjMask_("wallAdjMask", fvMsh_),
    wallDist_("wallDist", fvMsh_)
{
    if (solverDict.found("ImmersedBoundary"))
    {
        dictionary IBDict = solverDict.subDict("ImmersedBoundary");

        int nEntries = IBDict.size();

        xiStabilityFactor_ = IBDict.lookupOrDefault("xiStabilityFactor", 0.5);

        // Add shapes to IB according to dictionary entries
        for (int e = 0; e < nEntries; e++)
        {
            word entry = IBDict.toc()[e];

            if (IBDict.isDict(entry))
            {
                dictionary entryDict = IBDict.subDict(entry);

                if (word(entryDict.lookup("type")) == "cylinder")
                {
                    shapes_.append
                    (
                        new cylinder
                        (
                            vector(entryDict.lookup("start")),
                            vector(entryDict.lookup("end")),
                            readScalar(entryDict.lookup("radius")),
                            bool(entryDict.lookupOrDefault("inverted", false))
                        )
                    );
                }
                else if (word(entryDict.lookup("type")) == "sphere")
                {
                    shapes_.append
                    (
                        new sphere
                        (
                            vector(entryDict.lookup("center")),
                            readScalar(entryDict.lookup("radius")),
                            bool(entryDict.lookupOrDefault("inverted", false))
                        )
                    );
                }
                else
                {
                    FatalError
                        << "Unknown immersed boundary shape type "
                        << word(entryDict.lookup("type"))
                        << endl;
                    FatalError.exit();
                }
            }
        }

        // Cell centers
        const meshField<vector,MeshType>& CC =
            fvMsh_.metrics<MeshType>().cellCenters();

        // Set IB masks
        forAll(mask_, l)
        {
            forAll(mask_[l], d)
            {
                forAllCells(mask_[l][d], i, j, k)
                {
                    // Base block indices
                    const label x = i+1 - mask_[l][d].A().left();
                    const label y = j+1 - mask_[l][d].A().bottom();
                    const label z = k+1 - mask_[l][d].A().aft();

                    mask_[l][d](i,j,k) = 0.0;
                    wallAdjMask_[l][d](i,j,k) = 0.0;

                    if (this->isInside(CC[l][d](i,j,k)))
                    {
                        mask_[l][d](i,j,k) = 1.0;
                    }
                    else
                    {
                        vector c(CC[l][d](i,j,k));

                        if
                        (
                            this->isInside(CC[l][d].B()(x-1,y,z))
                        )
                        {
                            wallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(CC[l][d](i-1,j,k));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            wallDist_[l][d](i,j,k).left() = xi;
                        }
                        else
                        {
                            wallDist_[l][d](i,j,k).left() = -1.0;
                        }

                        if
                        (
                            this->isInside(CC[l][d].B()(x,y-1,z))
                        )
                        {
                            wallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(CC[l][d](i,j-1,k));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            wallDist_[l][d](i,j,k).bottom() = xi;
                        }
                        else
                        {
                            wallDist_[l][d](i,j,k).bottom() = -1.0;
                        }

                        if
                        (
                            this->isInside(CC[l][d].B()(x,y,z-1))
                        )
                        {
                            wallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(CC[l][d](i,j,k-1));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            wallDist_[l][d](i,j,k).aft() = xi;
                        }
                        else
                        {
                            wallDist_[l][d](i,j,k).aft() = -1.0;
                        }

                        if
                        (
                            this->isInside(CC[l][d].B()(x+1,y,z))
                        )
                        {
                            wallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(CC[l][d](i+1,j,k));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            wallDist_[l][d](i,j,k).right() = xi;
                        }
                        else
                        {
                            wallDist_[l][d](i,j,k).right() = -1.0;
                        }

                        if
                        (
                            this->isInside(CC[l][d].B()(x,y+1,z))
                        )
                        {
                            wallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(CC[l][d](i,j+1,k));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            wallDist_[l][d](i,j,k).top() = xi;
                        }
                        else
                        {
                            wallDist_[l][d](i,j,k).top() = -1.0;
                        }

                        if
                        (
                            this->isInside(CC[l][d].B()(x,y,z+1))
                        )
                        {
                            wallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(CC[l][d](i,j,k+1));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            wallDist_[l][d](i,j,k).fore() = xi;
                        }
                        else
                        {
                            wallDist_[l][d](i,j,k).fore() = -1.0;
                        }
                    }
                }
            }
        }
    }
    else
    {
        mask_ = Zero;
        wallAdjMask_ = Zero;
        wallDist_ = Zero;
    }
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
        if (dist == -1)
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

defineTemplateTypeNameAndDebug(immersedBoundary<colocated>, 0);
defineTemplateTypeNameAndDebug(immersedBoundary<staggered>, 0);

template class immersedBoundary<colocated>;
template class immersedBoundary<staggered>;

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam
