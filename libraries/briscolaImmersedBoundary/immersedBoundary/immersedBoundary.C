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

immersedBoundary::immersedBoundary
(
    IOdictionary solverDict,
    const fvMesh& fvMsh
)
:
    fvMsh_(fvMsh),
    colMask_("colMask", fvMsh_),
    stagMask_("stagMask", fvMsh_),
    colWallAdjMask_("colWallAdjMask", fvMsh_),
    stagWallAdjMask_("stagWallAdjMask", fvMsh_),
    colWallDist_("colWallDist", fvMsh_),
    stagWallDist_("stagWallDist", fvMsh_)
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

        // Colocated cell centers
        const colocatedVectorField& colCC =
            fvMsh_.metrics<colocated>().cellCenters();

        // Staggered cell centers
        const staggeredVectorField& stagCC =
            fvMsh_.metrics<staggered>().cellCenters();

        // Set colocated IB masks
        forAll(colMask_, l)
        {
            forAllCells(colMask_[l][0], i, j, k)
            {
                colMask_[l][0](i,j,k) = 0.0;
                colWallAdjMask_[l][0](i,j,k) = 0.0;

                if (this->isInside(colCC[l][0](i,j,k)))
                {
                    colMask_[l][0](i,j,k) = 1.0;
                }
                else
                {
                    vector c(colCC[l][0](i,j,k));

                    if
                    (
                        this->isInside(colCC[l][0](max(colMask_[l][0].A().left(), i-1),j,k))
                    )
                    {
                        colWallAdjMask_[l][0](i,j,k) = 1.0;

                        vector nb(colCC[l][0](i-1,j,k));
                        scalar wd = this->wallDistance(c, nb);
                        scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                        colWallDist_[l][0](i,j,k).left() = xi;
                    }
                    else
                    {
                        colWallDist_[l][0](i,j,k).left() = -1.0;
                    }

                    if
                    (
                        this->isInside(colCC[l][0](i,max(colMask_[l][0].A().bottom(), j-1),k))
                    )
                    {
                        colWallAdjMask_[l][0](i,j,k) = 1.0;

                        vector nb(colCC[l][0](i,j-1,k));
                        scalar wd = this->wallDistance(c, nb);
                        scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                        colWallDist_[l][0](i,j,k).bottom() = xi;
                    }
                    else
                    {
                        colWallDist_[l][0](i,j,k).bottom() = -1.0;
                    }

                    if
                    (
                        this->isInside(colCC[l][0](i,j,max(colMask_[l][0].A().aft(),k-1)))
                    )
                    {
                        colWallAdjMask_[l][0](i,j,k) = 1.0;

                        vector nb(colCC[l][0](i,j,k-1));
                        scalar wd = this->wallDistance(c, nb);
                        scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                        colWallDist_[l][0](i,j,k).aft() = xi;
                    }
                    else
                    {
                        colWallDist_[l][0](i,j,k).aft() = -1.0;
                    }

                    if
                    (
                        this->isInside(colCC[l][0](min(colMask_[l][0].A().right(), i+1),j,k))
                    )
                    {
                        colWallAdjMask_[l][0](i,j,k) = 1.0;

                        vector nb(colCC[l][0](i+1,j,k));
                        scalar wd = this->wallDistance(c, nb);
                        scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                        colWallDist_[l][0](i,j,k).right() = xi;
                    }
                    else
                    {
                        colWallDist_[l][0](i,j,k).right() = -1.0;
                    }

                    if
                    (
                        this->isInside(colCC[l][0](i,min(colMask_[l][0].A().top(), j+1),k))
                    )
                    {
                        colWallAdjMask_[l][0](i,j,k) = 1.0;

                        vector nb(colCC[l][0](i,j+1,k));
                        scalar wd = this->wallDistance(c, nb);
                        scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                        colWallDist_[l][0](i,j,k).top() = xi;
                    }
                    else
                    {
                        colWallDist_[l][0](i,j,k).top() = -1.0;
                    }

                    if
                    (
                        this->isInside(colCC[l][0](i,j,min(colMask_[l][0].A().fore(), k+1)))
                    )
                    {
                        colWallAdjMask_[l][0](i,j,k) = 1.0;

                        vector nb(colCC[l][0](i,j,k+1));
                        scalar wd = this->wallDistance(c, nb);
                        scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                        colWallDist_[l][0](i,j,k).fore() = xi;
                    }
                    else
                    {
                        colWallDist_[l][0](i,j,k).fore() = -1.0;
                    }
                }
            }
        }

        // Set staggered IB masks
        forAll(stagMask_, l)
        {
            forAll(stagMask_[l], d)
            {
                forAllCells(stagMask_[l][d], i, j, k)
                {
                    stagMask_[l][d](i,j,k) = 0.0;
                    stagWallAdjMask_[l][d](i,j,k) = 0.0;

                    if (this->isInside(stagCC[l][d](i,j,k)))
                    {
                        stagMask_[l][d](i,j,k) = 1.0;
                    }
                    else
                    {
                        vector c(stagCC[l][d](i,j,k));

                        if
                        (
                            this->isInside(stagCC[l][d](max(stagMask_[l][d].A().left(), i-1),j,k))
                        )
                        {
                            stagWallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(stagCC[l][d](i-1,j,k));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            stagWallDist_[l][d](i,j,k).left() = xi;
                        }
                        else
                        {
                            stagWallDist_[l][d](i,j,k).left() = -1.0;
                        }

                        if
                        (
                            this->isInside(stagCC[l][d](i,max(stagMask_[l][d].A().bottom(), j-1),k))
                        )
                        {
                            stagWallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(stagCC[l][d](i,j-1,k));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            stagWallDist_[l][d](i,j,k).bottom() = xi;
                        }
                        else
                        {
                            stagWallDist_[l][d](i,j,k).bottom() = -1.0;
                        }

                        if
                        (
                            this->isInside(stagCC[l][d](i,j,max(stagMask_[l][d].A().aft(),k-1)))
                        )
                        {
                            stagWallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(stagCC[l][d](i,j,k-1));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            stagWallDist_[l][d](i,j,k).aft() = xi;
                        }
                        else
                        {
                            stagWallDist_[l][d](i,j,k).aft() = -1.0;
                        }

                        if
                        (
                            this->isInside(stagCC[l][d](min(stagMask_[l][d].A().right(), i+1),j,k))
                        )
                        {
                            stagWallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(stagCC[l][d](i+1,j,k));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            stagWallDist_[l][d](i,j,k).right() = xi;
                        }
                        else
                        {
                            stagWallDist_[l][d](i,j,k).right() = -1.0;
                        }

                        if
                        (
                            this->isInside(stagCC[l][d](i,min(stagMask_[l][d].A().top(), j+1),k))
                        )
                        {
                            stagWallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(stagCC[l][d](i,j+1,k));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            stagWallDist_[l][d](i,j,k).top() = xi;
                        }
                        else
                        {
                            stagWallDist_[l][d](i,j,k).top() = -1.0;
                        }

                        if
                        (
                            this->isInside(stagCC[l][d](i,j,min(stagMask_[l][d].A().fore(), k+1)))
                        )
                        {
                            stagWallAdjMask_[l][d](i,j,k) = 1.0;

                            vector nb(stagCC[l][d](i,j,k+1));
                            scalar wd = this->wallDistance(c, nb);
                            scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                            stagWallDist_[l][d](i,j,k).fore() = xi;
                        }
                        else
                        {
                            stagWallDist_[l][d](i,j,k).fore() = -1.0;
                        }
                    }
                }
            }
        }
    }
    else
    {
        colMask_ = Zero;
        stagMask_ = Zero;
        colWallAdjMask_ = Zero;
        stagWallAdjMask_ = Zero;
        colWallDist_ = Zero;
        stagWallDist_ = Zero;
    }
}

// Destructor

immersedBoundary::~immersedBoundary()
{}

bool immersedBoundary::isInside(vector xyz)
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

scalar immersedBoundary::wallDistance(vector c, vector nb)
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

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam
