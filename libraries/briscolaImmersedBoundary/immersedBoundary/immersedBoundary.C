#include "immersedBoundary.H"

namespace Foam
{

namespace briscola
{

namespace ibm
{

using Foam::max;
using Foam::min;
using Foam::sqr;

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

        if (xiStabilityFactor_ < 0)
        {
            WarningInFunction
                << "Stability factor for immersed boundary"
                << " set to 0 (0 <= stability factor <= 0.999)"
                << endl;

            xiStabilityFactor_ = 0;
        }
        else if (xiStabilityFactor_ > 0.999)
        {
            WarningInFunction
                << "Stability factor for immersed boundary"
                << " set to 0.99 (0 <= stability factor <= 0.999)"
                << endl;

            xiStabilityFactor_ = 0.999;
        }

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
                    const label x = i+1 - mask_[l][d].I().left();
                    const label y = j+1 - mask_[l][d].I().bottom();
                    const label z = k+1 - mask_[l][d].I().aft();

                    mask_[l][d](i,j,k) = 0.0;
                    wallAdjMask_[l][d](i,j,k) = 0.0;
                    wallDist_[l][d](i,j,k) = -1.0;

                    if (this->isInside(CC[l][d](i,j,k)))
                    {
                        mask_[l][d](i,j,k) = 1.0;
                    }
                    else
                    {
                        vector c(CC[l][d](i,j,k));

                        for (int dir = 1; dir < 7; dir++)
                        {
                            // dir and faceOffsets are shifted by 1 due to
                            // stencil starting with center coefficient
                            const labelVector fo = faceOffsets[dir-1];

                            const labelVector xyzOffset =
                                labelVector(x,y,z) + fo;

                            if
                            (
                                this->isInside(CC[l][d].B()(xyzOffset))
                            )
                            {
                                wallAdjMask_[l][d](i,j,k) = 1.0;

                                vector nb(CC[l][d].B()(xyzOffset));
                                scalar wd = this->wallDistance(c, nb);
                                scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                                wallDist_[l][d](i,j,k)[dir] = xi;
                            }
                            else
                            {
                                wallDist_[l][d](i,j,k)[dir] = -1.0;
                            }
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

template<class MeshType>
void immersedBoundary<MeshType>::penalization
(
    linearSystem<stencil,scalar,MeshType>& ls
)
{
    // Set coefficients and sources to 0 in IB
    ls.A() *= (1.0 - mask_);
    ls.b() *= (1.0 - mask_);

    forAll(ls.A(), l)
    {
        forAll(ls.A()[l], d)
        {
            forAllCells(ls.A()[l][d], i, j, k)
            {
                // Set central coefficients to 1 in IB
                ls.A()[l][d](i,j,k).center() += mask_[l][d](i,j,k);
            }
        }
    }
}

template<class MeshType>
void immersedBoundary<MeshType>::IBM
(
    linearSystem<stencil,scalar,MeshType>& ls
)
{
    // We only apply the IBM on the finest mesh (l == 0)
    forAll(ls.A()[0], d)
    {
        forAllCells(ls.A()[0][d], i, j, k)
        {
            // Modify stencils in IB-adjacent cells
            if (wallAdjMask_[0][d](i,j,k) == 1)
            {
                // Loop over stencil directions (skipping center)
                for (int dir = 1; dir < 7; dir++)
                {
                    // dir and faceOffsets are shifted by 1 due to
                    // stencil starting with center coefficient
                    const label oppositeDir =
                        faceNumber(-faceOffsets[dir-1]) + 1;

                    if (wallDist_[0][d](i,j,k)[dir] >= 0)
                    {
                        scalar xi = wallDist_[0][d](i,j,k)[dir];

                        scalar w0 = 2.0 /
                            (
                                (1.0 - xiStabilityFactor_)
                                * (2.0 - xiStabilityFactor_)
                            );
                        scalar w1 = 2.0 - (2.0 - xi) * w0;
                        scalar w2 = -1.0 + (1.0 - xi) * w0;

                        if (xi < xiStabilityFactor_)
                        {
                            w1 = -2.0*xi/(1.0-xi);
                            w2 = xi/(2.0-xi);
                        }

                        const scalar a0 = ls.A()[0][d](i,j,k)[dir];

                        ls.A()[0][d](i,j,k)[dir] = 0;

                        ls.A()[0][d](i,j,k).center() += a0*w1;

                        ls.A()[0][d](i,j,k)[oppositeDir] += a0*w2;
                    }
                }
            }
        }
    }
}

template<class MeshType>
void immersedBoundary<MeshType>::imPCorr
(
    linearSystem<stencil,scalar,MeshType>& ls
)
{
    // We only apply the IBM on the finest mesh (l == 0)
    forAll(ls.A()[0], d)
    {
        forAllCells(ls.A()[0][d], i, j, k)
        {
            // Modify stencils in IB-adjacent cells
            if (wallAdjMask_[0][d](i,j,k) == 1)
            {
                // Loop over stencil directions (skipping center)
                for (int dir = 1; dir < 7; dir++)
                {
                    if (wallDist_[0][d](i,j,k)[dir] >= 0)
                    {
                        const scalar a = ls.A()[0][d](i,j,k)[dir];

                        ls.A()[0][d](i,j,k)[dir] = 0;

                        ls.A()[0][d](i,j,k).center() += a;
                    }
                }
            }
        }
    }
}

template<class MeshType>
tmp<colocatedScalarField> immersedBoundary<MeshType>::exPCorr
(
    colocatedScalarField& p
)
{
    tmp<colocatedScalarField> tpCorr
    (
        new colocatedScalarField
        (
            "tpCorr",
            p.fvMsh()
        )
    );

    colocatedScalarField& pCorr = tpCorr.ref();

    const PtrList<scalarList>& cellSizes
    (
        p.fvMsh().msh().cast<rectilinearMesh>().globalCellSizes()
    );

    forAllCells(pCorr[0][0], i, j, k)
    {
        pCorr[0][0](i,j,k) = 0;

        // Modify stencils in IB-adjacent cells
        if (wallAdjMask_[0][0](i,j,k) == 1)
        {
            if (wallDist_[0][0](i,j,k).left() >= 0)
            {
                pCorr[0][0](i,j,k) += p[0][0](i,j,k)
                    / sqr(cellSizes[0][i]);
                pCorr[0][0](i,j,k) -= p[0][0](i-1,j,k)
                    / sqr(cellSizes[0][i]);
            }

            if (wallDist_[0][0](i,j,k).bottom() >= 0)
            {
                pCorr[0][0](i,j,k) += p[0][0](i,j,k)
                    / sqr(cellSizes[1][j]);
                pCorr[0][0](i,j,k) -= p[0][0](i,j-1,k)
                    / sqr(cellSizes[1][j]);
            }

            if (wallDist_[0][0](i,j,k).aft() >= 0)
            {
                pCorr[0][0](i,j,k) += p[0][0](i,j,k)
                    / sqr(cellSizes[2][k]);
                pCorr[0][0](i,j,k) -= p[0][0](i,j,k-1)
                    / sqr(cellSizes[2][k]);
            }

            if (wallDist_[0][0](i,j,k).right() >= 0)
            {
                pCorr[0][0](i,j,k) += p[0][0](i,j,k)
                    / sqr(cellSizes[0][i]);
                pCorr[0][0](i,j,k) -= p[0][0](i+1,j,k)
                    / sqr(cellSizes[0][i]);
            }

            if (wallDist_[0][0](i,j,k).top() >= 0)
            {
                pCorr[0][0](i,j,k) += p[0][0](i,j,k)
                    / sqr(cellSizes[1][j]);
                pCorr[0][0](i,j,k) -= p[0][0](i,j+1,k)
                    / sqr(cellSizes[1][j]);
            }

            if (wallDist_[0][0](i,j,k).fore() >= 0)
            {
                pCorr[0][0](i,j,k) += p[0][0](i,j,k)
                    / sqr(cellSizes[2][k]);
                pCorr[0][0](i,j,k) -= p[0][0](i,j,k+1)
                    / sqr(cellSizes[2][k]);
            }
        }
    }

    return tpCorr;
}

defineTemplateTypeNameAndDebug(immersedBoundary<colocated>, 0);
defineTemplateTypeNameAndDebug(immersedBoundary<staggered>, 0);

template class immersedBoundary<colocated>;
template class immersedBoundary<staggered>;

} // end namespace ibm

} // end namespace briscola

} // end namespace Foam
