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
    type_
    (
        solverDict.found("ImmersedBoundary") ?
        solverDict.subDict("ImmersedBoundary").
            lookupOrDefault<word>("type", "Fadlun") : "none"
    ),
    massSourceActive_
    (
        solverDict.found("ImmersedBoundary") ?
        solverDict.subDict("ImmersedBoundary").
            lookupOrDefault<bool>("massSourceActive", false) : false
    ),
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
    ghostMask_
    (
        "ghostMask",
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
    // ,sendRecvList_(Pstream::nProcs())
{
    Info << "Immersed boundary type: " << type_ << endl;

    if (solverDict.found("ImmersedBoundary"))
    {
        // Set send/recv list sizes
        // for (int pr = 0; pr < Pstream::nProcs(); pr++)
        // {
        //     sendRecvList_[pr].setSize(Pstream::nProcs());

        //     for (int ps = 0; ps < Pstream::nProcs(); ps++)
        //     {
        //         const label nLevels = mask_.numberOfLevels();
        //         const label nDirections = MeshType::numberOfDirections;

        //         sendRecvList_[pr][ps].setSize(nLevels);

        //         for (int l = 0; l < nLevels; l++)
        //         {
        //             sendRecvList_[pr][ps][l].setSize(nDirections);

        //             for (int d = 0; d < nDirections; d++)
        //             {
        //                 sendRecvList_[pr][ps][l][d].setSize(0);
        //             }
        //         }
        //     }
        // }

        dictionary IBDict = solverDict.subDict("ImmersedBoundary");

        int nEntries = IBDict.size();

        // Add shapes to IB according to dictionary entries
        for (int e = 0; e < nEntries; e++)
        {
            word entry = IBDict.toc()[e];

            if (IBDict.isDict(entry))
            {
                dictionary entryDict = IBDict.subDict(entry);

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
        forAllLevels(mask_,l,d,i,j,k)
        {
            const labelVector ijk(i,j,k);

            mask_(l,d,i,j,k) = 0;
            wallAdjMask_(l,d,i,j,k) = 0;
            ghostMask_(l,d,i,j,k) = 0;
            wallDist_(l,d,i,j,k) = -1.0;
            neighborDist_(l,d,i,j,k) = -1.0;

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
                        wallAdjMask_(l,d,i,j,k) = 1;

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

        // Set recv lists
        // const decomposition& procDecomp = fvMsh_.msh().decomp();
        // forAllLevels(ghostMask_,l,d,i,j,k)
        // {
        //     if (ghostMask_(l,d,i,j,k) == 1)
        //     {
        //         const labelVector ijk(i,j,k);
        //         for (int dir = 0; dir < 6; dir++)
        //         {
        //             int xyzDir = dir / 2;
        //             const labelVector secondNeighbor
        //                 = ijk + 2.0*faceOffsets[dir];

        //             // If the neigboring cell in dir direction is a fluid cell
        //             // and the second neighbor is outside of the processor part
        //             if
        //             (
        //                 (mask_[l][d](ijk+faceOffsets[dir]) == 0)
        //                 &&
        //                 (
        //                     (dir%2 == 0) ?
        //                       (secondNeighbor[xyzDir] < ghostMask_[l][d].I()[dir])
        //                     : (secondNeighbor[xyzDir] > ghostMask_[l][d].I()[dir])
        //                 )
        //             )
        //             {
        //                 labelVector myProcIndex
        //                     = procDecomp.map().legend()[Pstream::myProcNo()];

        //                 labelVector neighborProcIndex = myProcIndex;
        //                 neighborProcIndex[xyzDir] +=  (dir%2 == 0) ? -1 : 1;

        //                 int neighborProc = procDecomp.map()(neighborProcIndex);

        //                 // Second neighbor global index
        //                 const labelVector globalSecondNeighbor = secondNeighbor
        //                     + procDecomp.globalStartPerProc()[Pstream::myProcNo()];

        //                 // Add index to list
        //                 sendRecvList_[Pstream::myProcNo()][neighborProc][l][d]
        //                     .append(globalSecondNeighbor);

        //                 break;
        //             }
        //         }
        //     }
        // }

        Pstream::gatherList(sendRecvList_);
        Pstream::scatterList(sendRecvList_);

    }
    else
    {
        mask_ = Zero;
        wallAdjMask_ = Zero;
        ghostMask_ = Zero;
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
scalar immersedBoundary<Type,MeshType>::wallNormalDistance
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

template<class Type, class MeshType>
vector immersedBoundary<Type,MeshType>::mirrorPoint
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

template<class Type, class MeshType>
void immersedBoundary<Type,MeshType>::penalization
(
    linearSystem<stencil,Type,MeshType>& ls
)
{
    forAllLevels(ls.b(),l,d,i,j,k)
    {
        if (mask_(l,d,i,j,k) == 1)
        {
            // Set sources to 0 in IB
            ls.b()(l,d,i,j,k) = Zero;
        }
    }

    forAllLevels(ls.A(),l,d,i,j,k)
    {
        if (mask_(l,d,i,j,k) == 1)
        {
            // Set coefficients to 0 in IB
            ls.A()(l,d,i,j,k) = Zero;

            // Set central coefficients to 1 in IB
            ls.A()(l,d,i,j,k).center() = 1.0;
        }
    }
}

template<class Type, class MeshType>
void immersedBoundary<Type,MeshType>::FadlunIBM
(
    linearSystem<stencil,Type,MeshType>& ls
)
{
    forAllLevels(ls.b(),l,d,i,j,k)
    {
        if (wallAdjMask_(l,d,i,j,k) == 1)
        {
            ls.b()(l,d,i,j,k) = Zero;
        }
    }

    forAllLevels(ls.A(),l,d,i,j,k)
    {
        // Modify stencils in IB-adjacent cells
        if (wallAdjMask_(l,d,i,j,k) == 1)
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

template<>
tmp<colocatedScalarField> immersedBoundary<scalar,staggered>::IBMSource
(
    const staggeredScalarField& field
)
{
    tmp<colocatedScalarField> tSource
    (
        new colocatedScalarField
        (
            "IBMSource",
            field.fvMsh()
        )
    );

    colocatedScalarField& source = tSource.ref();

    source = Zero;

    const colocatedFaceScalarField& fa =
        fvMsh_.metrics<colocated>().faceAreas();

    const colocatedScalarField& cv =
        fvMsh_.metrics<colocated>().cellVolumes();

    forAllCells(source[0][0],i,j,k)
    {
        source(0,0,i,j,k) -=
            ghostMask_(0,0,i,j,k) * field(0,0,i,j,k) * fa(0,0,i,j,k).left();

        source(0,0,i,j,k) +=
            ghostMask_(0,0,i+1,j,k) * field(0,0,i+1,j,k) * fa(0,0,i,j,k).right();

        source(0,0,i,j,k) -=
            ghostMask_(0,1,i,j,k) * field(0,1,i,j,k) * fa(0,0,i,j,k).bottom();

        source(0,0,i,j,k) +=
            ghostMask_(0,1,i,j+1,k) * field(0,1,i,j+1,k) * fa(0,0,i,j,k).top();

        source(0,0,i,j,k) -=
            ghostMask_(0,2,i,j,k) * field(0,2,i,j,k) * fa(0,0,i,j,k).aft();

        source(0,0,i,j,k) +=
            ghostMask_(0,2,i,j,k+1) * field(0,2,i,j,k+1) * fa(0,0,i,j,k).fore();

        source(0,0,i,j,k) /= cv(0,0,i,j,k);
    }

    return tSource;
}

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
