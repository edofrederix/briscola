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

template<class Type, class MeshType>
immersedBoundary<Type,MeshType>::immersedBoundary
(
    IOdictionary solverDict,
    const fvMesh& fvMsh
)
:
    fvMsh_(fvMsh),
    mask_
    (
        "mask",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        true
    ),
    wallAdjMask_
    (
        "wallAdjMask",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        true
    ),
    wallDist_
    (
        "wallDist",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        true
    ),
    neighborDist_
    (
        "neighborDist",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        false,
        true
    )
{
    if (solverDict.found("ImmersedBoundary"))
    {
        dictionary IBDict = solverDict.subDict("ImmersedBoundary");

        int nEntries = IBDict.size();

        xiStabilityFactor_
            = IBDict.lookupOrDefault("xiStabilityFactor", 0.999);

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
                << " set to 0.999 (0 <= stability factor <= 0.999)"
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

        // Set IB mask fields
        forAllLevels(mask_,l,d,i,j,k)
        {
            const labelVector ijk(i,j,k);

            mask_(l,d,i,j,k) = 0.0;
            wallAdjMask_(l,d,i,j,k) = 0.0;
            wallDist_(l,d,i,j,k) = -1.0;
            neighborDist_(l,d,i,j,k) = -1.0;

            if (this->isInside(CC(l,d,i,j,k)))
            {
                mask_(l,d,i,j,k) = 1.0;
            }
            else
            {
                const vector c(CC(l,d,i,j,k));

                for (int dir = 0; dir < 6; dir++)
                {
                    const labelVector fo = faceOffsets[dir];
                    const label oppositeDir = faceNumber(-fo);

                    if
                    (
                        this->isInside(CC[l][d](ijk+fo))
                    )
                    {
                        wallAdjMask_(l,d,i,j,k) = 1.0;

                        // Neighbor cell in the IB
                        const vector nb(CC[l][d](ijk+fo));
                        const scalar wd = this->wallDistance(c,nb);
                        const scalar xi = (mag(c-nb)-wd)/mag(c-nb);

                        wallDist_(l,d,i,j,k)[dir] = xi;

                        // Neighbor cell in the opposite direction
                        const vector nbo(CC[l][d](ijk-fo));
                        const scalar xi2 = mag(nbo-nb)/mag(c-nb);

                        neighborDist_(l,d,i,j,k)[oppositeDir] = xi2;
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
        neighborDist_ = Zero;
    }
}

template<class Type, class MeshType>
immersedBoundary<Type,MeshType>::~immersedBoundary()
{}

template<class Type, class MeshType>
bool immersedBoundary<Type,MeshType>::isInside(vector xyz)
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

template<class Type, class MeshType>
scalar immersedBoundary<Type,MeshType>::wallDistance(vector c, vector nb)
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

template<class Type, class MeshType>
void immersedBoundary<Type,MeshType>::penalization
(
    linearSystem<stencil,Type,MeshType>& ls
)
{
    // Set coefficients and sources to 0 in IB
    ls.A() *= (1.0 - mask_);
    ls.b() *= (1.0 - mask_);

    forAllLevels(ls.A(),l,d,i,j,k)
    {
        // Set central coefficients to 1 in IB
        ls.A()(l,d,i,j,k).center() += mask_(l,d,i,j,k);
    }
}

template<class Type, class MeshType>
void immersedBoundary<Type,MeshType>::IBM
(
    linearSystem<stencil,Type,MeshType>& ls
)
{
    forAllLevels(ls.A(),l,d,i,j,k)
    {
        // Modify stencils in IB-adjacent cells
        if (mag(wallAdjMask_(l,d,i,j,k) - 1.0) < 0.01)
        {
            // Loop over face number directions
            for (int dir = 0; dir < 6; dir++)
            {
                const label oppositeDir =
                    faceNumber(-faceOffsets[dir]);

                if (wallDist_(l,d,i,j,k)[dir] >= 0)
                {
                    const scalar xi
                        = wallDist_(l,d,i,j,k)[dir];

                    const scalar xi2
                        = neighborDist_(l,d,i,j,k)[oppositeDir];

                    scalar w0 = 2.0 /
                        (
                            (1.0 - xiStabilityFactor_)
                            * (2.0 - xiStabilityFactor_)
                        );
                    scalar w1 = 2.0 - (2.0 - xi) * w0;
                    scalar w2 = -1.0 + (1.0 - xi) * w0;

                    if (xi < xiStabilityFactor_)
                    {
                        w1 = xi*xi2/((1.0-xi)*(1.0-xi2));
                        w2 = xi/((xi2-xi)*(xi2-1.0));
                    }

                    // faceScalar and stencil directions are offset by 1

                    const scalar a0 = ls.A()(l,d,i,j,k)[dir+1];

                    ls.A()(l,d,i,j,k)[dir+1] = 0.0;

                    ls.A()(l,d,i,j,k).center() += a0*w1;

                    ls.A()(l,d,i,j,k)[oppositeDir+1] += a0*w2;
                }
            }
        }
    }
}

template<class Type, class MeshType>
void immersedBoundary<Type,MeshType>::IBM2
(
    linearSystem<stencil,Type,MeshType>& ls
)
{
    ls.b() *= (1.0 - wallAdjMask_);

    forAllLevels(ls.A(),l,d,i,j,k)
    {
        // Modify stencils in IB-adjacent cells
        if (mag(wallAdjMask_(l,d,i,j,k) - 1.0) < 0.01)
        {
            scalar ximax = 0;
            // Loop over face number directions
            for (int dir = 0; dir < 6; dir++)
            {
                const label oppositeDir =
                    faceNumber(-faceOffsets[dir]);

                if (wallDist_(l,d,i,j,k)[dir] > ximax)
                {
                    ximax = wallDist_(l,d,i,j,k)[dir];
                    const scalar xic = 1.0 - wallDist_(l,d,i,j,k)[dir];
                    const scalar xinb = 1.0 + xic;
                    const scalar w = xic/xinb;

                    ls.A()(l,d,i,j,k) = Zero;
                    ls.A()(l,d,i,j,k).center() = 1.0;
                    ls.A()(l,d,i,j,k)[oppositeDir+1] = -w;
                }
            }
        }
    }
}

typedef immersedBoundary<scalar,colocated> colocatedScalarImmersedBoundary;
typedef immersedBoundary<vector,colocated> colocatedVectorImmersedBoundary;
typedef immersedBoundary<scalar,staggered> staggeredScalarImmersedBoundary;
typedef immersedBoundary<vector,staggered> staggeredVectorImmersedBoundary;

defineTemplateTypeNameAndDebug(colocatedScalarImmersedBoundary, 0);
defineTemplateTypeNameAndDebug(colocatedVectorImmersedBoundary, 0);
defineTemplateTypeNameAndDebug(staggeredScalarImmersedBoundary, 0);
defineTemplateTypeNameAndDebug(staggeredVectorImmersedBoundary, 0);

template class immersedBoundary<scalar,colocated>;
template class immersedBoundary<vector,colocated>;
template class immersedBoundary<scalar,staggered>;
template class immersedBoundary<vector,staggered>;

} // end namespace fv

} // end namespace briscola

} // end namespace Foam
