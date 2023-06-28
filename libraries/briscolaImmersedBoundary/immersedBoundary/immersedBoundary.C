#include "immersedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace ibm
{

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
    stagWallAdjMask_("stagWallAdjMask", fvMsh_)
{
    if (solverDict.found("ImmersedBoundary"))
    {
        int nEntries = solverDict.subDict("ImmersedBoundary").size();

        // Add shapes to IB according to dictionary entries
        for (int e = 0; e < nEntries; e++)
        {
            word entry = solverDict.subDict("ImmersedBoundary").toc()[e];

            dictionary entryDict
                = solverDict.subDict("ImmersedBoundary").subDict(entry);

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

        // Colocated cell centers
        const colocatedVectorField& colCC =
            fvMsh_.metrics<colocated>().cellCenters();

        // Staggered cell centers
        const staggeredVectorField& stagCC =
            fvMsh_.metrics<staggered>().cellCenters();

        // Set colocated IB mask
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
                else if
                (
                        (i > colMask_[l][0].A().left() && this->isInside(colCC[l][0](i-1,j,k)))
                    || (j > colMask_[l][0].A().bottom() && this->isInside(colCC[l][0](i,j-1,k)))
                    || (k > colMask_[l][0].A().aft() && this->isInside(colCC[l][0](i,j,k-1)))
                    || (i < colMask_[l][0].A().right() - 1 && this->isInside(colCC[l][0](i+1,j,k)))
                    || (j < colMask_[l][0].A().top() - 1 && this->isInside(colCC[l][0](i,j+1,k)))
                    || (k < colMask_[l][0].A().fore() - 1 && this->isInside(colCC[l][0](i,j,k+1)))
                )
                {
                    colWallAdjMask_[l][0](i,j,k) = 1.0;
                }
            }
        }

        // Set staggered IB mask
        forAll(stagMask_, l)
        {
            forAll(stagMask_[l], d)
            {
                forAllCells(stagMask_[l][d], i, j, k)
                {
                    stagMask_[l][d](i,j,k) = 0.0;
                    stagMask_[l][d](i,j,k) = 0.0;

                    if (this->isInside(stagCC[l][d](i,j,k)))
                    {
                        stagMask_[l][d](i,j,k) = 1.0;
                    }
                    else
                    {
                        if
                        (
                               (i > stagMask_[l][0].A().left() && this->isInside(stagCC[l][d](i-1,j,k)))
                            || (j > stagMask_[l][0].A().bottom() && this->isInside(stagCC[l][d](i,j-1,k)))
                            || (k > stagMask_[l][0].A().aft() && this->isInside(stagCC[l][d](i,j,k-1)))
                            || (i < stagMask_[l][0].A().right() - 1 && this->isInside(stagCC[l][d](i+1,j,k)))
                            || (j < stagMask_[l][0].A().top() - 1 && this->isInside(stagCC[l][d](i,j+1,k)))
                            || (k < stagMask_[l][0].A().fore() - 1 && this->isInside(stagCC[l][d](i,j,k+1)))
                        )
                        {
                            stagWallAdjMask_[l][d](i,j,k) = 1.0;
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

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam
