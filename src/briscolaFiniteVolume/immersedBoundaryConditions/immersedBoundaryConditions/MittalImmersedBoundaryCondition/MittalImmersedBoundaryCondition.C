#include "MittalImmersedBoundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

template<class Type, class MeshType>
MittalImmersedBoundaryCondition<Type,MeshType>::MittalImmersedBoundaryCondition
(
    const meshField<Type,MeshType>& field,
    const immersedBoundary<MeshType>& ib
)
:
    immersedBoundaryCondition<Type,MeshType>(field, ib, &ib.ghostMask()),
    exchanges_(field.fvMsh().msh().size())
{
    // Check shape overlap

    if (this->ib_.shapeOverlap())
        this->shapeOverlapWarning();

    // Set mirror points

    const meshField<vector,MeshType>& mps = this->ib_.mirrorPoints();
    const meshField<label,MeshType>& mask = this->forcingMask();

    forAll(mask, l)
    {
        exchanges_.set
        (
            l,
            new PtrList<pointDataExchange<MeshType>>
            (
                MeshType::numberOfDirections
            )
        );

        forAll(mask[l], d)
        {
            label nPoints = 0;

            forAllCells(mask[l][d],i,j,k)
                if (mask(l,d,i,j,k))
                    nPoints++;

            vectorList exchangePoints(nPoints);

            label c = 0;

            forAllCells(mask[l][d],i,j,k)
            if (mask(l,d,i,j,k))
            {
                vector mp = mps(l,d,i,j,k);

                if (this->ib_.isInside(mp))
                {
                    WarningInFunction
                        << "Mirror point of ghost cell " << vector(i,j,k)
                        << " at level " << l << ", direction " << d
                        << " located inside of immersed boundary." << endl;
                }

                exchangePoints[c++] = mp;
            }

            exchanges_[l].set
            (
                d,
                new pointDataExchange<MeshType>
                (
                    exchangePoints,
                    this->fvMsh_,
                    l,
                    d
                )
            );
        }
    }
}

// Destructor

template<class Type, class MeshType>
MittalImmersedBoundaryCondition<Type,MeshType>::
~MittalImmersedBoundaryCondition()
{}

}

}

}
